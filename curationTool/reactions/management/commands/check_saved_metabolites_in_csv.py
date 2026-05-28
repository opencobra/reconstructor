"""
Check which saved metabolites (per user) can be found in the unmatched VMH CSV
by matching InChI keys at a configurable layer depth.
"""
import os
import sys
import django
from collections import defaultdict

# ---------------------------------------------------------------------------
# Globals
# ---------------------------------------------------------------------------
CSV_PATH = "/home/saleh/Downloads/2026_04_20_unmapped_metabolites_v7.xlsx - Unmatched VMH.csv"

# InChI key layer depth (number of dash-separated segments to compare):
#   1 = connectivity layer only  (e.g. "UXFQDXABPXWSTK")
#   2 = + stereo layer           (e.g. "UXFQDXABPXWSTK-ZDUSSCGKSA")
#   3 = full key (default)       (e.g. "UXFQDXABPXWSTK-ZDUSSCGKSA-L")
INCHI_LEVEL = 3


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def truncate_inchikey(key: str, level: int) -> str:
    """Return the first `level` dash-separated segments of an InChI key."""
    if not key:
        return ""
    parts = key.strip().split("-")
    return "-".join(parts[:level])


def load_csv_inchikeys(csv_path: str, level: int) -> set:
    """Return the set of (truncated) InChI keys present in the CSV."""
    import csv

    keys: set = set()
    with open(csv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            raw = row.get("InChIKey", "").strip()
            if raw:
                keys.add(truncate_inchikey(raw, level))
    return keys


def setup_django():
    repo_root = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    )
    sys.path.insert(0, repo_root)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "reconstructor.settings")
    django.setup()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    setup_django()

    from reactions.models import SavedMetabolite

    print(f"Loading CSV: {CSV_PATH}")
    print(f"InChI key match level: {INCHI_LEVEL} segment(s)\n")

    csv_keys = load_csv_inchikeys(CSV_PATH, INCHI_LEVEL)
    print(f"  {len(csv_keys)} unique InChI keys loaded from CSV.\n")

    # Group saved metabolites by user
    user_found: dict[str, list] = defaultdict(list)    # user_label -> [(name, id, inchikey)]
    user_missing: dict[str, list] = defaultdict(list)  # user_label -> [(name, id, inchikey)]

    qs = SavedMetabolite.objects.select_related("owner").only(
        "id", "name", "inchi_key", "owner__id", "owner__username"
    )

    for met in qs.iterator():
        owner = met.owner
        user_label = (
            f"{owner.username} (id={owner.id})" if owner else "Unknown"
        )
        raw_key = met.inchi_key or ""
        truncated = truncate_inchikey(raw_key, INCHI_LEVEL)

        entry = (met.name or "", met.id, raw_key)

        if truncated and truncated in csv_keys:
            user_found[user_label].append(entry)
        else:
            user_missing[user_label].append(entry)

    # ---------------------------------------------------------------------------
    # Report
    # ---------------------------------------------------------------------------
    all_users = sorted(set(user_found) | set(user_missing))
    total_found = 0
    total_missing = 0

    print("=" * 70)
    print("Per-user summary")
    print("=" * 70)

    for user_label in all_users:
        found = user_found.get(user_label, [])
        missing = user_missing.get(user_label, [])
        total_found += len(found)
        total_missing += len(missing)

        print(f"\n{user_label}")
        print(f"  Found in CSV : {len(found)}")
        print(f"  NOT in CSV   : {len(missing)}")

        if missing:
            print("  --- Not found ---")
            for name, met_id, inchi_key in missing:
                print(f"    id={met_id:<6}  [{inchi_key or 'no key'}]  {name}")

    print("\n" + "=" * 70)
    print(f"TOTAL  found={total_found}  not_found={total_missing}")
    print("=" * 70)


if __name__ == "__main__":
    main()
