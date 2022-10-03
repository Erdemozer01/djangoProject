# Generated by Django 4.1.1 on 2022-10-03 09:44

import ckeditor_uploader.fields
from django.db import migrations, models
import django.utils.timezone


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='About',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(default='Hakkımızda', max_length=100, verbose_name='Hakkımızda')),
                ('about', ckeditor_uploader.fields.RichTextUploadingField(verbose_name='Açıklama:')),
                ('publish', models.DateTimeField(default=django.utils.timezone.now, verbose_name='Yayınlama Tarihi')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('updated', models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')),
            ],
            options={
                'verbose_name': 'Hakkımızda',
                'verbose_name_plural': 'Hakkımızda',
                'db_table': 'about',
            },
        ),
        migrations.CreateModel(
            name='Bottom',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(help_text='Anasayfa', max_length=100, verbose_name='Başlık')),
                ('text', models.TextField(verbose_name='İçerik')),
                ('publish', models.DateTimeField(default=django.utils.timezone.now, verbose_name='Yayınlama Tarihi')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('updated', models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')),
            ],
            options={
                'verbose_name': 'Alt Kapak',
                'verbose_name_plural': 'Alt Kapak',
                'db_table': 'bottom_cover',
                'ordering': ['publish'],
            },
        ),
        migrations.CreateModel(
            name='Contact',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(max_length=20, verbose_name='Başlık')),
                ('email', models.EmailField(max_length=254, verbose_name='Email')),
                ('phone', models.CharField(max_length=13, verbose_name='Telefon')),
                ('address', models.TextField(verbose_name='Adres')),
                ('publish', models.DateTimeField(default=django.utils.timezone.now, verbose_name='Yayınlama Tarihi')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('updated', models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')),
            ],
            options={
                'verbose_name': 'İletişim',
                'verbose_name_plural': 'İletişim',
                'db_table': 'contact',
                'ordering': ['publish'],
            },
        ),
        migrations.CreateModel(
            name='Cover',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('image', models.ImageField(blank=True, upload_to='home/cover/', verbose_name='Anasayfa Kapak')),
                ('title', models.CharField(help_text='Anasayfa', max_length=100, verbose_name='Başlık')),
                ('text', models.TextField(verbose_name='İçerik')),
                ('theme', models.CharField(blank=True, choices=[('', '------'), ('dark', 'Siyah'), ('primary', 'Pembe'), ('danger', 'Kırmızı'), ('info', 'Mavi'), ('success', 'Yaşil')], default='', help_text='Kapak Fotoğrafı seçtiyseniz seçin', max_length=100, verbose_name='Tema')),
                ('publish', models.DateTimeField(default=django.utils.timezone.now, verbose_name='Yayınlama Tarihi')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('updated', models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')),
            ],
            options={
                'verbose_name': 'Üst Kapak',
                'verbose_name_plural': 'Üst Kapak',
                'db_table': 'cover',
                'ordering': ['publish'],
            },
        ),
        migrations.CreateModel(
            name='Inbox',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=20, verbose_name='Ad Soyad:')),
                ('email', models.EmailField(max_length=254, verbose_name='Email:')),
                ('topic', models.CharField(max_length=50, verbose_name='Konu:')),
                ('phone', models.CharField(max_length=13, verbose_name='Telefon')),
                ('message', models.TextField(verbose_name='Mesaj')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Gönderilme Tarihi')),
            ],
            options={
                'verbose_name': 'Gelen Mesajlar',
                'verbose_name_plural': 'Gelen Mesajlar',
                'db_table': 'inbox',
                'ordering': ['created'],
            },
        ),
        migrations.CreateModel(
            name='Social',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(help_text='ör: Facebook', max_length=20, verbose_name='Sosyal Medya')),
                ('url', models.URLField(verbose_name='URL')),
                ('publish', models.DateTimeField(default=django.utils.timezone.now, verbose_name='Yayınlama Tarihi')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('updated', models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')),
            ],
            options={
                'verbose_name': 'Sosyal Medya',
                'verbose_name_plural': 'Sosyal Medya',
                'db_table': 'social',
            },
        ),
        migrations.CreateModel(
            name='Terms',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=100, verbose_name='Koşul ve Şartlar')),
                ('terms', ckeditor_uploader.fields.RichTextUploadingField(verbose_name='Koşul ve Şartlar')),
                ('publish', models.DateTimeField(default=django.utils.timezone.now, verbose_name='Gönderilme Tarihi')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('updated', models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')),
            ],
            options={
                'verbose_name': 'Koşul ve Şartlar',
                'verbose_name_plural': 'Koşul ve Şartlar',
                'ordering': ['created'],
            },
        ),
        migrations.CreateModel(
            name='Titles',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=20, verbose_name='Site Adı')),
                ('type', models.CharField(max_length=20, verbose_name='Site Türü')),
                ('publish', models.DateTimeField(default=django.utils.timezone.now, verbose_name='Yayınlama Tarihi')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('updated', models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')),
            ],
            options={
                'verbose_name': 'Site Başlığı',
                'verbose_name_plural': 'Site Başlığı',
                'db_table': 'title',
            },
        ),
    ]
