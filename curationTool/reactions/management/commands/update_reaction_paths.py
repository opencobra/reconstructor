from django.core.management.base import BaseCommand
from django.core.exceptions import ObjectDoesNotExist
import json
from reactions.models import Reaction  # Adjust the import path as necessary

class Command(BaseCommand):
    help = 'Updates the visualization field in Reaction objects to remove "reactions/" from paths'

    def handle(self, *args, **options):
        updated_count = 0
        for reaction in Reaction.objects.all():
            try:
                visualization_data = json.loads(reaction.visualization)
                # Check if visualization is a list and not empty
                if isinstance(visualization_data, list) and visualization_data:
                    original_path = visualization_data[0]
                    # Update the first path if it contains 'reactions/'
                    if original_path.startswith('reactions/'):
                        new_path = original_path.replace('reactions/', '', 1)
                        visualization_data[0] = new_path
                        reaction.visualization = json.dumps(visualization_data)
                        reaction.save()
                        updated_count += 1
                        # Updated line to output detailed update information
                        self.stdout.write(self.style.SUCCESS(f'Updated Reaction ID {reaction.id}: path from "{original_path}" to "{new_path}"'))
            except json.JSONDecodeError:
                self.stdout.write(self.style.WARNING(f'Skipping Reaction ID {reaction.id} due to invalid JSON in visualization field'))
            except Exception as e:
                self.stdout.write(self.style.ERROR(f'Error updating Reaction ID {reaction.id}: {str(e)}'))

        self.stdout.write(self.style.SUCCESS(f'Finished updating {updated_count} reactions.'))
