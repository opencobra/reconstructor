import os
import django
from django.core.management.base import BaseCommand
from reactions.models import User, CreatedReaction

class Command(BaseCommand):
    help = 'Copy saved reactions to CreatedReaction for each user.'

    def handle(self, *args, **kwargs):
        self.copy_saved_reactions_to_created_reactions()

    def copy_saved_reactions_to_created_reactions(self):
        users = User.objects.all()
        for user in users:
            saved_reactions = user.saved_reactions.all()
            for reaction in saved_reactions:
                # Check if the CreatedReaction already exists to avoid duplication
                created_reaction, created = CreatedReaction.objects.get_or_create(
                    user=user,
                    reaction=reaction,
                )
                if created:
                    self.stdout.write(self.style.SUCCESS(f"CreatedReaction for user {user.name} and reaction {reaction.short_name} created."))
                else:
                    self.stdout.write(self.style.WARNING(f"CreatedReaction for user {user.name} and reaction {reaction.short_name} already exists."))
