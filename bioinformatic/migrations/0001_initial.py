# Generated by Django 4.1.1 on 2022-09-24 00:03

import bioinformatic.models
from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='BigFileUploadModel',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('big_file', models.FileField(upload_to='laboratory/bigfile/', verbose_name='Dosya Seçme')),
                ('created', models.DateTimeField(auto_now_add=True, null=True)),
            ],
            options={
                'verbose_name': 'bigfile',
                'verbose_name_plural': 'bigfile',
                'db_table': 'bigfile',
            },
        ),
        migrations.CreateModel(
            name='BlastQuery',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('blast', models.TextField()),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
            ],
            options={
                'ordering': ['-created'],
            },
        ),
        migrations.CreateModel(
            name='FastaRead',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('gene', models.CharField(max_length=1000)),
                ('sequence', models.TextField()),
            ],
            options={
                'verbose_name': 'Fasta Okuması',
                'verbose_name_plural': 'Fasta Okumaları',
                'db_table': 'fasta_read',
            },
        ),
        migrations.CreateModel(
            name='FileUploadModel',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('file', models.FileField(upload_to='laboratory/bigfile/', verbose_name='Dosya Seçme')),
                ('created', models.DateTimeField(auto_now_add=True, null=True)),
            ],
            options={
                'verbose_name': 'file',
                'verbose_name_plural': 'files',
                'db_table': 'file',
            },
        ),
        migrations.CreateModel(
            name='GenbankRead',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('organism', models.CharField(max_length=1000, verbose_name='Organizma')),
                ('description', models.CharField(max_length=1000, verbose_name='Tanım')),
                ('protein_id', models.CharField(blank=True, max_length=255, null=True, verbose_name='protein_id')),
                ('taxonomy', models.CharField(blank=True, max_length=1000, null=True, verbose_name='Taksonomi')),
                ('dna_sequence', models.TextField()),
                ('dna_sequence_len', models.BigIntegerField(blank=True, null=True, verbose_name='DNA Uzunluğu')),
                ('protein_sequence', models.TextField(blank=True, null=True)),
                ('protein_sequence_len', models.BigIntegerField(blank=True, null=True, verbose_name='Protein Uzunluğu')),
            ],
            options={
                'verbose_name': 'Genbank Okuması',
                'verbose_name_plural': 'Genbank Okumaları',
                'db_table': 'genbank_read',
                'ordering': ['id'],
            },
        ),
        migrations.CreateModel(
            name='LabSlideModel',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('image', models.ImageField(upload_to='laboratory/slide/')),
                ('title', models.CharField(max_length=100)),
                ('text', models.CharField(max_length=100)),
                ('content', models.TextField()),
            ],
            options={
                'verbose_name': 'Slide',
                'verbose_name_plural': 'Slide',
                'db_table': 'bioinformatic_home',
            },
        ),
        migrations.CreateModel(
            name='MedlineArticle',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('email', models.EmailField(max_length=254)),
                ('article_id', models.CharField(max_length=1000)),
                ('doi', models.CharField(default='', max_length=1000)),
                ('title', models.CharField(max_length=1000)),
                ('term', models.CharField(max_length=1000)),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
            ],
        ),
        migrations.CreateModel(
            name='PubMedArticle',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('email', models.EmailField(max_length=254)),
                ('article_id', models.CharField(max_length=100)),
                ('title', models.CharField(max_length=1000)),
                ('link', models.URLField(default='https://pubmed.ncbi.nlm.nih.gov/')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
            ],
            options={
                'verbose_name': 'PubMed Makale',
                'verbose_name_plural': 'PubMed Makale',
            },
        ),
        migrations.CreateModel(
            name='SwissProtModel',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('accessions', models.CharField(max_length=1000, verbose_name='Erişim Numarası')),
                ('organism', models.CharField(max_length=1000, verbose_name='Organizma')),
                ('sequence', models.TextField(verbose_name='Sekans')),
                ('sequence_length', models.CharField(max_length=1000, verbose_name='Sekans Uzunluğu')),
            ],
            options={
                'verbose_name': 'swiss-prot',
                'verbose_name_plural': 'swiss-prot',
                'db_table': 'swiss-prot',
            },
        ),
        migrations.CreateModel(
            name='MultipleSequenceAlignment',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('method', models.CharField(choices=[('', '------------'), ('MUSCLE', 'MUSCLE'), ('clustalw2', 'CLUSTALW2'), ('omega', 'ClustalOmega')], max_length=1000, verbose_name='Multiple Sekans Alignment Metodu Seçiniz')),
                ('algoritma', models.CharField(choices=[('', '------------'), ('nj', 'Neighbor Joining'), ('upgma', 'UPGMA')], max_length=1000, verbose_name='Algoritma Seçiniz')),
                ('molecule_type', models.CharField(choices=[('', '------------'), ('DNA', 'DNA'), ('RNA', 'RNA'), ('protein', 'PROTEİN')], max_length=1000, verbose_name='Molekül Tipi')),
                ('alignment_filetype', models.CharField(choices=[('', '------------'), ('fasta', 'FASTA'), ('clustal', 'CLUSTAL'), ('phylib', 'PHYLİB'), ('nexus', 'NEXUS')], max_length=1000, verbose_name='Alignment Dosya Tipi')),
                ('input_file', models.FileField(upload_to=bioinformatic.models.upload_to, verbose_name='Fasta Dosyası')),
                ('output_file', models.FileField(blank=True, null=True, upload_to=bioinformatic.models.upload_to, verbose_name='Alignment Dosyası')),
                ('align_file', models.FileField(blank=True, null=True, upload_to=bioinformatic.models.upload_to, verbose_name='Aligned Dosyası')),
                ('stats', models.FileField(blank=True, null=True, upload_to=bioinformatic.models.upload_to, verbose_name='İstatistik Dosyası')),
                ('scores', models.FileField(blank=True, null=True, upload_to=bioinformatic.models.upload_to, verbose_name='Alignment Scor Dosyası')),
                ('tree', models.ImageField(upload_to='', verbose_name='Filogenetik Ağaç')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL, verbose_name='Laborant')),
            ],
            options={
                'verbose_name': 'Multiple Sequence Alignment',
                'verbose_name_plural': 'Multiple Sequence Alignment',
                'db_table': 'msa',
            },
        ),
    ]
