# Generated by Django 4.1.3 on 2022-11-30 02:02

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bioinformatic', '0008_molecularmodel_author_molecularmodel_head_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='molecularmodel',
            name='file_name',
            field=models.CharField(blank=True, max_length=100, null=True, verbose_name='Dosya Adı:'),
        ),
    ]