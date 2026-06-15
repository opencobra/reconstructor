from django.core.management.base import BaseCommand
from reactions.models import Subsystem


class Command(BaseCommand):
    help = "Load subsystem names into the local Subsystem table (deduplicated against existing rows)"

    def add_arguments(self, parser):
        parser.add_argument(
            "filepath",
            type=str,
            help="Path to a text file with one subsystem name per line",
        )

    def handle(self, *args, **kwargs):
        filepath = kwargs["filepath"]

        try:
            with open(filepath, encoding="utf-8") as f:
                names = [line.strip() for line in f]
        except FileNotFoundError:
            self.stderr.write(self.style.ERROR(f"File not found: {filepath}"))
            return

        created = 0
        skipped = 0
        for name in names:
            if not name:  # skip blank / whitespace-only lines (e.g. the empty VMH row)
                continue
            _, was_created = Subsystem.objects.get_or_create(name=name)
            if was_created:
                created += 1
            else:
                skipped += 1

        self.stdout.write(self.style.SUCCESS(
            f"Done. {created} new subsystem(s) added, {skipped} already existed."
        ))
