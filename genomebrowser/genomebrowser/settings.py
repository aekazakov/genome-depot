"""
Django settings for genomebrowser project.

Generated by 'django-admin startproject' using Django 3.0.5.

For more information on this file, see
https://docs.djangoproject.com/en/3.0/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/3.0/ref/settings/
"""

import os
#import json
from decouple import config, Csv

# Normally you should not import ANYTHING from Django directly
# into your settings, but ImproperlyConfigured is an exception.
from django.core.exceptions import ImproperlyConfigured

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# JSON-based secrets module
#with open(os.path.join(BASE_DIR, "secrets.json")) as f:
#    secrets = json.loads(f.read())

#def get_secret(setting, secrets=secrets):
#    """Get the secret variable or return explicit exception."""
#    try:
#        return secrets[setting]
#    except KeyError:
#        error_msg = "Set the {0} environment variable".format(setting)
#        raise ImproperlyConfigured(error_msg)


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/3.0/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production in the .env file!
SECRET_KEY = config('SECRET_KEY')

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = config('DEBUG', cast=bool)

ALLOWED_HOSTS = config('ALLOWED_HOSTS', cast=Csv())
INTERNAL_IPS = config('INTERNAL_IPS', cast=Csv()) #['127.0.0.1',]
TITLE = config('TITLE')

# Application definition

INSTALLED_APPS = [
    'browser.apps.BrowserConfig',
    'admin_shortcuts',
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_q',
    'debug_toolbar',
]

if DEBUG:
    INSTALLED_APPS += ('corsheaders', )

CORS_ORIGIN_ALLOW_ALL = DEBUG

MIDDLEWARE = [
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'debug_toolbar.middleware.DebugToolbarMiddleware',
]

ROOT_URLCONF = 'genomebrowser.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'genomebrowser.wsgi.application'


# Database
# https://docs.djangoproject.com/en/3.0/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': config('DB_NAME'),
        'USER': config('DB_USER'),
        'PASSWORD': config('DB_PASSWORD'),
        'CONN_MAX_AGE': 3600,
        'TEST': {
            'NAME':'test_genomesdev'
        },
        'PORT': '3306',
        'HOST': '127.0.0.1',
        'OPTIONS': {
           'init_command': 'SET default_storage_engine=INNODB',
           'isolation_level': 'read committed'
        }
    }
}


# Password validation
# https://docs.djangoproject.com/en/3.0/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
{'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',},
{'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',},
{'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',},
{'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',},
]


# Internationalization
# https://docs.djangoproject.com/en/3.0/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'America/Los_Angeles'

USE_I18N = True

USE_L10N = True

USE_TZ = True

ADMINS = [tuple(config('ADMIN', cast=Csv()))]
EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'
EMAIL_USE_TLS = True
EMAIL_PORT = 587
EMAIL_HOST = 'smtp.gmail.com'
#sender's email-id
EMAIL_HOST_USER = config('EMAIL_HOST_USER')
#password associated with the email-id
EMAIL_HOST_PASSWORD = config('EMAIL_HOST_PASSWORD')

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/3.0/howto/static-files/

BASE_URL = config('BASE_URL')
STATIC_URL = config('STATIC_URL')
STATIC_ROOT = config('STATIC_ROOT')
STATICFILES_DIRS = [
    config('STATICFILES_DIR')
]

DEFAULT_AUTO_FIELD='django.db.models.AutoField'

CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'cgcms_cache_table',
    }
}

Q_CLUSTER = {
    'name': 'cgcmsWorker',
    'workers': 1,
    'timeout': 8640000,
    'retry': 8640001,
    'ack_failures': True,
    'save_limit': 0,
    'max_attempts': 1,
    'attempt_count': 1,
    'orm': 'default',
}

DEBUG_TOOLBAR_PANELS = [
    'debug_toolbar.panels.versions.VersionsPanel',
    'debug_toolbar.panels.timer.TimerPanel',
    'debug_toolbar.panels.settings.SettingsPanel',
    'debug_toolbar.panels.headers.HeadersPanel',
    'debug_toolbar.panels.request.RequestPanel',
    'debug_toolbar.panels.sql.SQLPanel',
    'debug_toolbar.panels.staticfiles.StaticFilesPanel',
    'debug_toolbar.panels.templates.TemplatesPanel',
    'debug_toolbar.panels.cache.CachePanel',
    'debug_toolbar.panels.signals.SignalsPanel',
    'debug_toolbar.panels.logging.LoggingPanel',
    'debug_toolbar.panels.redirects.RedirectsPanel',
]

ADMIN_SHORTCUTS = [
    {
        'title': 'Genomes and annotations',
        'shortcuts': [
            {
                'title': 'Genomes',
                'icon': 'folder',
                'url_name': 'admin:browser_genome_changelist',
            },
            {
                'title': 'Import genomes',
                'url_name': 'admin:browser_genome_add',
                'count': 'from genbank files',
            },
            {
                'title': 'Annotations',
                'icon': 'folder',
                'url_name': 'admin:browser_annotation_changelist',
            },
            {
                'title': 'Import annotations',
                'url_name': 'admin:browser_annotation_add',
                'url': 'browser/annotation/import-annotations',
                'count': 'from tab-separated file',
            },
            {
                'title': 'Regulons',
                'icon': 'folder',
                'url_name': 'admin:browser_regulon_changelist',
            },
            {
                'title': 'Import regulons',
                'url_name': 'admin:browser_regulon_add',
                'url': 'browser/regulon/import-regulons',
                'count': 'from tab-separated file',
            },
        ]
    },
    {
        'title': 'Strains and samples',
        'shortcuts': [
            {
                'title': 'Strains',
                'icon': 'folder',
                'url_name': 'admin:browser_strain_changelist',
            },
            {
                'title': 'Strain metadata',
                'icon': 'folder',
                'url_name': 'admin:browser_strain_metadata_changelist',
            },
            {
                'title': 'Update strain metadata',
                'url_name': 'admin:browser_strain_metadata_add',
                'count': 'from Excel file',
            },
            {
                'title': 'Samples',
                'icon': 'folder',
                'url_name': 'admin:browser_sample_changelist',
            },
            {
                'title': 'Import sample descriptions',
                'url_name': 'admin:browser_sample_add',
                'url': 'browser/sample/import-descriptions',
                'count': 'from tab-separated file',
            },
            {
                'title': 'Sample metadata',
                'icon': 'folder',
                'url_name': 'admin:browser_sample_metadata_changelist',
            },
            {
                'title': 'Import sample metadata',
                'url_name': 'admin:browser_sample_metadata_add',
                'count': 'from tab-separated file',
            },
        ]
    },
    {
        'title': 'Tools',
        'shortcuts': [
            {
                'title': 'Tasks',
                'icon': 'truck',
                'url_name': 'admin:django_q_ormq_changelist',
                'count_new': 'browser.admin.count_tasks',
                'count': 'Check queue',
            },
            {
                'title': 'Configuration',
                'icon': 'edit',
                'url_name': 'admin:browser_config_changelist',
                'count': 'Edit configs',
            },
            {
                'title': 'Groups',
                'url_name': 'admin:auth_group_changelist',
                'icon': 'columns',
                'count': 'View user groups',
            },
            {
                'title': 'Add user',
                'url_name': 'admin:auth_user_add',
                'has_perms': 'example.utils.has_perms_to_users',
            },
            {
                'title': 'Clusters',
                'url': 'clusters/',
                'count_new': 'browser.admin.count_clusters',
                'count':'Django Q clusters',
            },
        ]
    },
]
