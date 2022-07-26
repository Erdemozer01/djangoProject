# Generated by Django 4.1.3 on 2022-11-27 15:27

import bioinformatic.models
from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('bioinformatic', '0003_alter_graphicmodels_format'),
    ]

    operations = [
        migrations.AlterField(
            model_name='multiplesequencealignment',
            name='method',
            field=models.CharField(choices=[('', '------------'), ('muscle', 'Muscle'), ('clustalw2', 'ClustalW2'), ('omega', 'Clustal Omega'), ('paml', 'Maximum Likelihood (PAML)')], max_length=1000, verbose_name='Multiple Sekans Alignment Aracı'),
        ),
        migrations.AlterField(
            model_name='multiplesequencealignment',
            name='palm_tools',
            field=models.CharField(blank=True, choices=[('', '------------'), ('baseml', 'BASEML'), ('basemlg', 'BASEMLG'), ('codeml', 'CODEML')], max_length=1000, null=True, verbose_name='Maximum Likelihood (PAML)'),
        ),
        migrations.CreateModel(
            name='MolecularModel',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('in_file', models.FileField(blank=True, null=True, upload_to=bioinformatic.models.molecule_upload_to, verbose_name='PDB Dosyası')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL, verbose_name='Laborant')),
            ],
        ),
    ]
