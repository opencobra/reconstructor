from django.core.management.base import BaseCommand
from django.contrib.auth.models import User

class Command(BaseCommand):
    help = 'Displays all admin accounts'

    def handle(self, *args, **kwargs):
        admin_users = User.objects.filter(is_superuser=True) | User.objects.filter(is_staff=True)
        for user in admin_users:
            self.stdout.write(f"Username: {user.username}, Superuser: {user.is_superuser}, Staff: {user.is_staff}, Password: {user.password}")
