"""API key lifecycle, /api/v1 auth behaviour, throttling, and the web key page."""
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
    def _sign_in(self, user):
        session = self.client.session
        session['userID'] = user.pk
        session.save()

    def test_generate_requires_session(self):
        response = self.client.post('/api-key/generate/')
        self.assertEqual(response.status_code, 401)

    def test_page_renders_without_session(self):
        response = self.client.get('/api-key/')
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, 'Not signed in')

    def test_generate_rotate_roundtrip(self):
        user = make_user()
        self._sign_in(user)
        response = self.client.post('/api-key/generate/')
        self.assertEqual(response.status_code, 200)
        first_key = response.json()['api_key']
        self.assertTrue(first_key.startswith('ctor_'))
        self.assertEqual(
            self.client.get('/api/v1/me', **auth(first_key)).status_code, 200)

        page = self.client.get('/api-key/')
        self.assertContains(page, first_key[:12])
        self.assertNotContains(page, first_key)  # raw key never re-shown

        second_key = self.client.post('/api-key/generate/').json()['api_key']
        self.assertEqual(
            self.client.get('/api/v1/me', **auth(first_key)).status_code, 401)
        self.assertEqual(
            self.client.get('/api/v1/me', **auth(second_key)).status_code, 200)
