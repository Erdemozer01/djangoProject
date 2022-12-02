# Generated by Django 4.1.3 on 2022-11-30 01:32

import bioinformatic.models
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bioinformatic', '0007_molecularmodel_file_path_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='molecularmodel',
            name='author',
            field=models.CharField(blank=True, max_length=100, null=True, verbose_name='Yazar:'),
        ),
        migrations.AddField(
            model_name='molecularmodel',
            name='head',
            field=models.CharField(blank=True, max_length=100, null=True, verbose_name='Özellik:'),
        ),
        migrations.AddField(
            model_name='molecularmodel',
            name='id_code',
            field=models.CharField(blank=True, max_length=100, null=True, verbose_name='id cod:'),
        ),
        migrations.AddField(
            model_name='molecularmodel',
            name='keywords',
            field=models.CharField(blank=True, max_length=100, null=True, verbose_name='Kategori:'),
        ),
        migrations.AddField(
            model_name='molecularmodel',
            name='name',
            field=models.CharField(blank=True, max_length=100, null=True, verbose_name='Adı:'),
        ),
        migrations.AlterField(
            model_name='molecularmodel',
            name='file_path',
            field=models.FilePathField(blank=True, null=True, path=bioinformatic.models.molecule_path, recursive=True),
        ),
        migrations.AlterField(
            model_name='molecularmodel',
            name='id_name',
            field=models.CharField(blank=True, max_length=100, null=True, verbose_name='id:'),
        ),
        migrations.AlterField(
            model_name='molecularmodel',
            name='in_file',
            field=models.FileField(blank=True, null=True, upload_to=bioinformatic.models.molecule_upload_to2, verbose_name='PDB Dosyası'),
        ),
    ]
