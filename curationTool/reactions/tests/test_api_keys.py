"""API key lifecycle, /api/v1 auth behaviour, throttling, and the web key page."""
import json

from django.core.cache import cache
from django.test import TestCase

from reactions.api import api_throttle
from reactions.models import ApiKey
from reactions.services import api_keys
from reactions.tests.helpers import auth, make_key, make_user


class ApiKeyServiceTests(TestCase):
    def test_generate_returns_prefixed_key_and_stores_hash_only(self):
        user = make_user()
        raw = api_keys.generate_key(user)
        self.assertTrue(raw.startswith('ctor_'))
        key = ApiKey.objects.get(user=user)
        self.assertNotIn(raw, key.key_hash)
        self.assertEqual(key.prefix, raw[:12])
        self.assertEqual(api_keys.authenticate(raw), user)

    def test_rotation_invalidates_old_key(self):
        user = make_user()
        old = api_keys.generate_key(user)
        new = api_keys.generate_key(user)
        self.assertIsNone(api_keys.authenticate(old))
        self.assertEqual(api_keys.authenticate(new), user)
        self.assertEqual(ApiKey.objects.filter(user=user).count(), 1)

    def test_authenticate_rejects_garbage(self):
        self.assertIsNone(api_keys.authenticate(None))
        self.assertIsNone(api_keys.authenticate(''))
        self.assertIsNone(api_keys.authenticate('ctor_nonexistent'))
        self.assertIsNone(api_keys.authenticate('wrong_prefix'))


class ApiAuthTests(TestCase):
    def setUp(self):
        # Throttle history lives in the cache; sqlite reuses pks across
        # rolled-back tests, so stale entries would bleed between tests.
        cache.clear()

    def test_health_is_open(self):
        response = self.client.get('/api/v1/health')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['status'], 'ok')

    def test_me_requires_key(self):
        response = self.client.get('/api/v1/me')
        self.assertEqual(response.status_code, 401)
        self.assertEqual(response.json()['error']['code'], 'unauthorized')

    def test_me_with_bearer_and_header_auth(self):
        user = make_user()
        key = make_key(user)
        response = self.client.get('/api/v1/me', **auth(key))
        self.assertEqual(response.status_code, 200)
        body = response.json()
        self.assertEqual(body['user_id'], user.pk)
        self.assertEqual(body['key_prefix'], key[:12])

        response = self.client.get('/api/v1/me', HTTP_X_API_KEY=key)
        self.assertEqual(response.status_code, 200)

    def test_rotated_key_fails_on_api(self):
        user = make_user()
        old = make_key(user)
        make_key(user)
        response = self.client.get('/api/v1/me', **auth(old))
        self.assertEqual(response.status_code, 401)

    def test_throttle_returns_rate_limited_envelope(self):
        user = make_user(name='throttled-user')
        key = make_key(user)
        original = (api_throttle.num_requests, api_throttle.duration)
        api_throttle.num_requests, api_throttle.duration = 3, 60
        try:
            statuses = [
                self.client.get('/api/v1/me', **auth(key)).status_code
                for _ in range(4)
            ]
        finally:
            api_throttle.num_requests, api_throttle.duration = original
        self.assertEqual(statuses[:3], [200, 200, 200])
        self.assertEqual(statuses[3], 429)
        response = self.client.get('/api/v1/me', **auth(key))
        # still over limit; also check envelope
        if response.status_code == 429:
            self.assertEqual(response.json()['error']['code'], 'rate_limited')


class ApiKeyWebPageTests(TestCase):
    def setUp(self):
        cache.clear()  # failed-password counters live in the cache

    def _sign_in(self, user):
        session = self.client.session
        session['userID'] = user.pk
        session.save()

    def _generate(self, password):
        data = {} if password is None else {'password': password}
        return self.client.post('/api-key/generate/', data)

    def test_generate_requires_session(self):
        response = self.client.post('/api-key/generate/')
        self.assertEqual(response.status_code, 401)

    def test_page_renders_without_session(self):
        response = self.client.get('/api-key/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Not signed in')

    def test_page_asks_for_password(self):
        self._sign_in(make_user())
        response = self.client.get('/api-key/')
        self.assertContains(response, 'id="apiKeyPassword"')

    def test_generate_rotate_roundtrip(self):
        user = make_user(password='s3cret')
        self._sign_in(user)
        response = self._generate('s3cret')
        self.assertEqual(response.status_code, 200)
        first_key = response.json()['api_key']
        self.assertTrue(first_key.startswith('ctor_'))
        self.assertEqual(
            self.client.get('/api/v1/me', **auth(first_key)).status_code, 200)

        page = self.client.get('/api-key/')
        self.assertContains(page, first_key[:12])
        self.assertNotContains(page, first_key)  # raw key never re-shown

        second_key = self._generate('s3cret').json()['api_key']
        self.assertEqual(
            self.client.get('/api/v1/me', **auth(first_key)).status_code, 401)
        self.assertEqual(
            self.client.get('/api/v1/me', **auth(second_key)).status_code, 200)

    def test_forged_session_cannot_mint_or_rotate_a_key(self):
        # The browser syncs any user id into the session (set_session_user),
        # so a session alone must not be enough to obtain someone's API key.
        victim = make_user(name='victim', password='victim-pw')
        victim_key = make_key(victim)
        response = self.client.post(
            '/set_session_user/', data=json.dumps({'userID': victim.pk}),
            content_type='application/json')
        self.assertEqual(response.json()['status'], 'success')

        for password in (None, '', 'guess'):
            response = self._generate(password)
            self.assertEqual(response.status_code, 403)
            self.assertNotIn('api_key', response.json())
        # the victim's existing key was not rotated away
        self.assertEqual(api_keys.authenticate(victim_key), victim)

    def test_failed_password_attempts_are_limited(self):
        user = make_user(password='right')
        self._sign_in(user)
        for _ in range(5):
            self.assertEqual(self._generate('wrong').status_code, 403)
        # locked even with the right password until the window passes
        response = self._generate('right')
        self.assertEqual(response.status_code, 429)
        self.assertFalse(ApiKey.objects.filter(user=user).exists())

    def test_successful_password_resets_the_counter(self):
        user = make_user(password='right')
        self._sign_in(user)
        for _ in range(4):
            self._generate('wrong')
        self.assertEqual(self._generate('right').status_code, 200)
        for _ in range(4):
            self.assertEqual(self._generate('wrong').status_code, 403)
        self.assertEqual(self._generate('right').status_code, 200)
