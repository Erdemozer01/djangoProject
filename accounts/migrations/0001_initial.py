# Generated by Django 4.1 on 2022-08-13 20:02

import accounts.models
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
            name='Profile',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('first_name', models.CharField(max_length=150, verbose_name='Ad')),
                ('last_name', models.CharField(max_length=150, verbose_name='Soyad')),
                ('email', models.EmailField(max_length=254, verbose_name='Email')),
                ('cover', models.ImageField(blank=True, upload_to=accounts.models.cover, verbose_name='Kapak Fotosu:')),
                ('avatar', models.ImageField(blank=True, upload_to=accounts.models.avatar, verbose_name='Avatar:')),
                ('gender', models.CharField(blank=True, choices=[('', '------'), ('Erkek', 'Erkek'), ('Kadın', 'Kadın'), ('Belirtmek istemiyorum', 'Belirtmek istemiyorum')], default='', max_length=30, verbose_name='Cinsiyet')),
                ('location', models.CharField(blank=True, max_length=100, verbose_name='Yaşadığı şehir')),
                ('about', models.TextField(blank=True, verbose_name='Hakkımda')),
                ('job', models.CharField(blank=True, max_length=100, verbose_name='Meslek')),
                ('phone', models.CharField(blank=True, max_length=13, null=True, unique=True, verbose_name='Telefon')),
                ('birth_day', models.DateField(blank=True, null=True, verbose_name='Doğum Tarihi:')),
                ('skils', models.CharField(blank=True, max_length=200, verbose_name='Yetenekler')),
                ('facebook', models.URLField(blank=True, verbose_name='Facebook')),
                ('twitter', models.URLField(blank=True, verbose_name='Twitter')),
                ('instagram', models.URLField(blank=True, verbose_name='İnstagram')),
                ('created', models.DateTimeField(auto_now=True, verbose_name='Katılma Tarihi')),
                ('user', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='profile', to=settings.AUTH_USER_MODEL, verbose_name='Kullanıcı Adı')),
            ],
            options={
                'verbose_name': 'Profil',
                'verbose_name_plural': 'Profil',
                'db_table': 'profile',
            },
        ),
    ]
