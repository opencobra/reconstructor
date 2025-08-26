import csv
import os
from django.core.management.base import BaseCommand
from reactions.models import Gene

class Command(BaseCommand):
    help = "Load HGNC gene dataset into database"

    def add_arguments(self, parser):
        parser.add_argument("filepath", type=str, help="Path to hgnc_complete_set.txt")

    def handle(self, *args, **kwargs):
        filepath = kwargs["filepath"]

        if not os.path.exists(filepath):
            self.stderr.write(self.style.ERROR(f"File not found: {filepath}"))
            return

        with open(filepath, newline='', encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            count = 0
            for row in reader:
                symbol = row.get("symbol")
                hgnc_id = row.get("hgnc_id")
                name = row.get("name", "")
                entrez_id = row.get("entrez_id") or None
                aliases = ",".join(row.get("alias_symbol", "").split("|")) if row.get("alias_symbol") else ""

                if not symbol or not hgnc_id:
                    continue

                Gene.objects.update_or_create(
                    hgnc_id=hgnc_id,
                    defaults={
                        "symbol": symbol,
                        "name": name,
                        "entrez_id": entrez_id,
                        "aliases": aliases,
                    },
                )
                count += 1

        self.stdout.write(self.style.SUCCESS(f"Loaded/updated {count} genes"))
