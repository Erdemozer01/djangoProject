# Generated by Django 4.1 on 2022-08-22 03:59

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bioinformatic', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='pubmedarticle',
            name='link',
            field=models.URLField(default=''),
        ),
    ]
