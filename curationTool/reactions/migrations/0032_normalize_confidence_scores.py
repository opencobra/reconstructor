from django.db import migrations


VALID_CONFIDENCE_SCORES = {'1', '2', '3', '4'}


def _normalize(raw):
    """Mirror of reactions.utils.utils.normalize_confidence_score (kept inline
    so the migration is self-contained and stable over time)."""
    if raw is None:
        return None
    s = str(raw).strip()
    while len(s) >= 2 and s[0] in '"\'' and s[-1] == s[0]:
        s = s[1:-1].strip()
    return s if s in VALID_CONFIDENCE_SCORES else None


def normalize_confidence_scores(apps, schema_editor):
    Reaction = apps.get_model('reactions', 'Reaction')
    for reaction in Reaction.objects.exclude(confidence_score__isnull=True).iterator():
        normalized = _normalize(reaction.confidence_score)
        if normalized != reaction.confidence_score:
            reaction.confidence_score = normalized
            reaction.save(update_fields=['confidence_score'])


def noop(apps, schema_editor):
    # Normalization is not reversible (original encodings are intentionally lost).
    pass


class Migration(migrations.Migration):

    dependencies = [
        ('reactions', '0031_reaction_vmh_added_via_constructor'),
    ]

    operations = [
        migrations.RunPython(normalize_confidence_scores, noop),
    ]
