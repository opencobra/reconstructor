# Generated by Django 5.0.1 on 2024-08-23 10:44

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Reaction',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('substrates', models.TextField(help_text='Comma-separated list of substrates.')),
                ('products', models.TextField(help_text='Comma-separated list of products.')),
                ('short_name', models.TextField(blank=True, null=True)),
                ('direction', models.TextField(blank=True, null=True)),
                ('substrates_types', models.TextField(blank=True, null=True)),
                ('products_types', models.TextField(blank=True, null=True)),
                ('substrates_names', models.TextField(blank=True, null=True)),
                ('products_names', models.TextField(blank=True, null=True)),
                ('formulas', models.TextField(blank=True, null=True)),
                ('visualization', models.TextField(blank=True, null=True)),
                ('molc_formula', models.TextField(blank=True, null=True)),
                ('balanced_count', models.TextField(blank=True, null=True)),
                ('balanced_charge', models.TextField(blank=True, null=True)),
                ('subsystem', models.TextField(blank=True, null=True)),
                ('subs_comps', models.TextField(blank=True, null=True)),
                ('prods_comps', models.TextField(blank=True, null=True)),
                ('subs_sch', models.TextField(blank=True, null=True)),
                ('prods_sch', models.TextField(blank=True, null=True)),
                ('subs_atoms', models.TextField(blank=True, null=True)),
                ('prods_atoms', models.TextField(blank=True, null=True)),
                ('subs_charge', models.TextField(blank=True, null=True)),
                ('prods_charge', models.TextField(blank=True, null=True)),
                ('symb_to_name', models.TextField(blank=True, null=True)),
                ('metabolite_names', models.TextField(blank=True, null=True)),
                ('metabolite_formulas', models.TextField(blank=True, null=True)),
                ('metabolite_charges', models.TextField(blank=True, null=True)),
                ('metabolite_mol_file_strings', models.TextField(blank=True, null=True)),
                ('subs_found', models.TextField(blank=True, null=True)),
                ('subs_miriams', models.TextField(blank=True, null=True)),
                ('prod_found', models.TextField(blank=True, null=True)),
                ('prod_miriams', models.TextField(blank=True, null=True)),
                ('subs_edited', models.TextField(blank=True, null=True)),
                ('prods_edited', models.TextField(blank=True, null=True)),
                ('vmh_found', models.BooleanField(default=False)),
                ('vmh_found_similar', models.BooleanField(default=False)),
                ('vmh_url', models.TextField(blank=True, null=True)),
                ('vmh_formula', models.TextField(blank=True, null=True)),
                ('references', models.JSONField(blank=True, null=True)),
                ('ext_links', models.JSONField(blank=True, null=True)),
                ('gene_info', models.JSONField(blank=True, null=True)),
                ('comments', models.JSONField(blank=True, null=True)),
                ('confidence_score', models.CharField(blank=True, max_length=10, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Subsystem',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255, unique=True)),
            ],
        ),
        migrations.CreateModel(
            name='User',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=100, null=True)),
                ('password', models.CharField(blank=True, max_length=128, null=True)),
                ('cred_add_to_vmh', models.BooleanField(default=False)),
                ('cred_add_to_rhea', models.BooleanField(default=False)),
                ('orchid_id', models.CharField(blank=True, max_length=255, null=True)),
                ('email', models.EmailField(blank=True, max_length=254, null=True)),
                ('full_name', models.CharField(blank=True, max_length=255, null=True)),
                ('saved_reactions', models.ManyToManyField(blank=True, related_name='saved_by_users', to='reactions.reaction')),
            ],
        ),
        migrations.CreateModel(
            name='ReactionsAddedVMH',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('user_name', models.CharField(max_length=255)),
                ('reaction_id', models.CharField(max_length=255)),
                ('reaction_formula', models.TextField()),
                ('reaction_abbr', models.CharField(max_length=255)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='reactions.user')),
            ],
        ),
        migrations.CreateModel(
            name='MetabolitesAddedVMH',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('user_name', models.CharField(max_length=255)),
                ('metabolite_id', models.CharField(max_length=255)),
                ('metabolite_formula', models.TextField()),
                ('metabolite_abbr', models.CharField(max_length=255)),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='reactions.user')),
            ],
        ),
        migrations.CreateModel(
            name='CreatedReaction',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(auto_now_add=True)),
                ('reaction', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='created_by_users', to='reactions.reaction')),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='created_reactions', to='reactions.user')),
            ],
        ),
    ]
