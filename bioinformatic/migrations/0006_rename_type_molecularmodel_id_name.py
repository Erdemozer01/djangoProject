# Generated by Django 4.1.3 on 2022-11-28 15:11

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('bioinformatic', '0005_molecularmodel_type'),
    ]

    operations = [
        migrations.RenameField(
            model_name='molecularmodel',
            old_name='type',
            new_name='id_name',
        ),
    ]