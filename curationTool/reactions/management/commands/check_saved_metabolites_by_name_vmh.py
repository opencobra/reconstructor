"""
Check which saved metabolites (per user) exist in VMH by name (exact match only).
"""
from collections import defaultdict

import urllib3
from django.core.management.base import BaseCommand
from tqdm import tqdm

from reactions.models import SavedMetabolite
from reactions.utils.vmh_api import search_metabolites_new

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


class Command(BaseCommand):
    help = (
        "Check which SavedMetabolites exist in VMH by exact name match, "
        "grouped per user."
    )

    @staticmethod
    def _exact_name_match(rows: list, name: str) -> bool:
        target = (name or "").strip()
        if not target:
            return False
        for row in rows:
            full_name = (row.get("fullName") or "").strip()
            if full_name == target:
                return True
        return False

    def handle(self, *args, **options):
        user_found: dict = defaultdict(list)    # user_label -> [(name, id)]
        user_missing: dict = defaultdict(list)  # user_label -> [(name, id)]

        qs = SavedMetabolite.objects.select_related("owner").only(
            "id", "name", "owner__id", "owner__name"
        )
        total = qs.count()

        for met in tqdm(qs.iterator(), total=total, desc="Checking VMH by name"):
            owner = met.owner
            user_label = f"{owner.name} (id={owner.id})" if owner else "Unknown"
            name = (met.name or "").strip()
            entry = (name, met.id)

            if not name:
                user_missing[user_label].append(entry)
                continue

            try:
                rows = search_metabolites_new({"fullName": name})
            except Exception as exc:
                self.stderr.write(f"Error searching '{name}' (id={met.id}): {exc}")
                user_missing[user_label].append(entry)
                continue

            if self._exact_name_match(rows, name):
                user_found[user_label].append(entry)
            else:
                user_missing[user_label].append(entry)

        all_users = sorted(set(user_found) | set(user_missing))
        total_found = 0
        total_missing = 0

        self.stdout.write("\n" + "=" * 70)
        self.stdout.write("Per-user summary (VMH exact name match)")
        self.stdout.write("=" * 70)

        for user_label in all_users:
            found = user_found.get(user_label, [])
            missing = user_missing.get(user_label, [])
            total_found += len(found)
            total_missing += len(missing)

            self.stdout.write(f"\n{user_label}")
            self.stdout.write(f"  Found in VMH : {len(found)}")
            self.stdout.write(f"  NOT in VMH   : {len(missing)}")

            if found:
                self.stdout.write("  --- Found ---")
                for name, met_id in found:
                    self.stdout.write(f"    id={met_id:<6}  {name}")
            if missing:
                self.stdout.write("  --- Not found ---")
                for name, met_id in missing:
                    self.stdout.write(f"    id={met_id:<6}  {name}")

        self.stdout.write("\n" + "=" * 70)
        self.stdout.write(f"TOTAL  found={total_found}  not_found={total_missing}")
        self.stdout.write("=" * 70)
