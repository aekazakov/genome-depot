from django.contrib.auth.models import User, Group

def has_perms_to_users(request):
    return request.user.is_authenticated and request.user.is_staff
