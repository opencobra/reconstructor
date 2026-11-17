"""Settings for running tests and local smoke servers without Postgres.

Usage:
    python manage.py test reactions --settings=reactions_project.test_settings
    python manage.py runserver --settings=reactions_project.test_settings
"""
import sys

from .settings import *  # noqa: F401,F403

TESTING = sys.argv[1:2] == ['test']

if TESTING:
    # Tests never talk to VMH: no key, an unresolvable API host (so an unmocked
    # call fails fast instead of reaching the real service), per-process caches
    # and no derived InChIKey index unless a test opts in.
    VMH_API_KEY = ''
    VMH_API_BASE_URL = 'http://vmh.invalid/api/public'
    VMH_BASE_URL = 'https://www.vmh.life'
    VMH_DERIVED_INDEX = False
    CHEBI_API_BASE_URL = 'http://chebi.invalid/chebi/backend/api/public'
    CACHES = {
        'default': {'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
                    'LOCATION': 'constructor-tests-default'},
        'vmh': {'BACKEND': 'django.core.cache.backends.locmem.LocMemCache',
                'LOCATION': 'constructor-tests-vmh'},
    }

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'test_db.sqlite3',  # noqa: F405
    }
}

# No collectstatic manifest in tests / local smoke runs.
STATICFILES_STORAGE = 'django.contrib.staticfiles.storage.StaticFilesStorage'

# Serve static/media from runserver so the sandbox web UI renders with its
# CSS/JS. Harmless for tests: Django's test runner forces DEBUG=False.
DEBUG = True
SERVE_MEDIA = True
