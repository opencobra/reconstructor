"""Settings for running tests and local smoke servers without Postgres.

Usage:
    python manage.py test reactions --settings=reactions_project.test_settings
    python manage.py runserver --settings=reactions_project.test_settings
"""
from .settings import *  # noqa: F401,F403

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': BASE_DIR / 'test_db.sqlite3',  # noqa: F405
    }
}

# No collectstatic manifest in tests / local smoke runs.
STATICFILES_STORAGE = 'django.contrib.staticfiles.storage.StaticFilesStorage'
