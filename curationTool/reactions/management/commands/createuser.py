from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError

class Command(BaseCommand):
    help = 'Creates a new user'

    def add_arguments(self, parser):
        parser.add_argument('username', type=str, help='Username for the new user')
        parser.add_argument('password', type=str, help='Password for the new user')
        parser.add_argument('--superuser', action='store_true', help='Create a superuser account')

    def handle(self, *args, **options):
        username = options['username']
        password = options['password']
        is_superuser = options['superuser']

        try:
            if is_superuser:
                User.objects.create_superuser(username=username, password=password)
            else:
                User.objects.create_user(username=username, password=password)
            self.stdout.write(self.style.SUCCESS(f'Successfully created user: {username}'))
        except ValidationError as e:
            self.stdout.write(self.style.ERROR(f'Error: {e}'))
