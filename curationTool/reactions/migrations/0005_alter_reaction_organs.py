# Generated by Django 5.0.1 on 2024-08-24 16:26

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('reactions', '0004_reaction_organs'),
    ]

    operations = [
        migrations.AlterField(
            model_name='reaction',
            name='Organs',
            field=models.JSONField(blank=True, null=True),
        ),
    ]