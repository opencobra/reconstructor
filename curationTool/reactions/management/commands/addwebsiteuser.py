from django.core.management.base import BaseCommand
from django.core.exceptions import ValidationError
from reactions.models import User  # Import your User model

"""python manage.py addwebsiteuser username password"""

class Command(BaseCommand):
    help = 'Creates a new website user ensuring unique key and username'

    def add_arguments(self, parser):
        parser.add_argument('--username', type=str, help='Username for the new user')
        parser.add_argument('--password', type=str, help='Password for the new user')

    def create_user(self, username, password):
        if User.objects.filter(name=username).exists():
            raise ValidationError("Username already exists.")

        User.objects.create(name=username, password=password)
        self.stdout.write(self.style.SUCCESS(f'Successfully created user: {username} with password: {password}'))

    def handle(self, *args, **options):
        username = options['username']
        password = options['password']
        try:
            self.create_user(username, password)
        except ValidationError as e:
            self.stdout.write(self.style.ERROR(f'Error: {e}'))
