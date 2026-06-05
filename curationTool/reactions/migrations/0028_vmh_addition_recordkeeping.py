# Generated manually for VMH addition record keeping.

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('reactions', '0027_gene_alter_savedmetabolite_inchi_key_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='added_metabolites',
            field=models.JSONField(blank=True, default=list),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='comment_snapshot',
            field=models.JSONField(blank=True, default=list),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='confidence_score',
            field=models.CharField(blank=True, default='', max_length=10),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='direction',
            field=models.CharField(blank=True, default='', max_length=64),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='external_link_snapshot',
            field=models.JSONField(blank=True, default=list),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='gene_info_snapshot',
            field=models.JSONField(blank=True, default=list),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='local_reaction',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='vmh_addition_records', to='reactions.reaction'),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='matlab_result',
            field=models.JSONField(blank=True, default=dict),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='product_snapshot',
            field=models.JSONField(blank=True, default=list),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='reaction_name',
            field=models.TextField(blank=True, default=''),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='reference_snapshot',
            field=models.JSONField(blank=True, default=list),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='request_snapshot',
            field=models.JSONField(blank=True, default=dict),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='substrate_snapshot',
            field=models.JSONField(blank=True, default=list),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='subsystem',
            field=models.TextField(blank=True, default=''),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='user_email',
            field=models.EmailField(blank=True, default='', max_length=254),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='user_full_name',
            field=models.CharField(blank=True, default='', max_length=255),
        ),
        migrations.AddField(
            model_name='reactionsaddedvmh',
            name='vmh_database',
            field=models.CharField(blank=True, default='', max_length=255),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='compartment',
            field=models.CharField(blank=True, default='', max_length=32),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='inchi_key',
            field=models.CharField(blank=True, default='', max_length=255),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='local_reaction',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='vmh_metabolite_addition_records', to='reactions.reaction'),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='metabolite_name',
            field=models.TextField(blank=True, default=''),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='reaction_abbr',
            field=models.CharField(blank=True, default='', max_length=255),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='reaction_entry',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='metabolite_records', to='reactions.reactionsaddedvmh'),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='reaction_formula',
            field=models.TextField(blank=True, default=''),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='side',
            field=models.CharField(blank=True, default='', max_length=32),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='source_identifier',
            field=models.TextField(blank=True, default=''),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='source_type',
            field=models.CharField(blank=True, default='', max_length=64),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='stoichiometry',
            field=models.CharField(blank=True, default='', max_length=64),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='submission_snapshot',
            field=models.JSONField(blank=True, default=dict),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='user_email',
            field=models.EmailField(blank=True, default='', max_length=254),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='user_full_name',
            field=models.CharField(blank=True, default='', max_length=255),
        ),
        migrations.AddField(
            model_name='metabolitesaddedvmh',
            name='vmh_database',
            field=models.CharField(blank=True, default='', max_length=255),
        ),
        migrations.AlterField(
            model_name='metabolitesaddedvmh',
            name='metabolite_formula',
            field=models.TextField(blank=True, default=''),
        ),
    ]
