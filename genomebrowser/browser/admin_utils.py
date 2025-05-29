from django.conf import settings

def has_perms_to_users(request):
    return request.user.is_authenticated and request.user.is_staff

def is_email_configured():
    required_settings = [
        'EMAIL_BACKEND',
        'EMAIL_HOST',
        'EMAIL_PORT',
        'EMAIL_HOST_USER',
        'EMAIL_HOST_PASSWORD',
    ]
    for setting in required_settings:
        if not hasattr(settings, setting) or not getattr(settings, setting):
            return False
    return True
