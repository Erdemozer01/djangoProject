# Generated by Django 4.1.1 on 2022-09-24 00:03

import autoslug.fields
import ckeditor_uploader.fields
from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Category',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('image', models.ImageField(upload_to='blog/category/', verbose_name='Kategori Fotosu:')),
                ('title', models.CharField(max_length=100, verbose_name='Kategori:')),
                ('explain', models.TextField(verbose_name='Kategori Tanımı:')),
                ('slug', autoslug.fields.AutoSlugField(editable=False, populate_from='title', unique=True)),
                ('publish', models.DateTimeField(default=django.utils.timezone.now, verbose_name='Yayınlama Tarihi')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('updated', models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')),
            ],
            options={
                'verbose_name': 'Kategori',
                'verbose_name_plural': 'Kategori',
                'db_table': 'category',
            },
        ),
        migrations.CreateModel(
            name='Posts',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('cover', models.ImageField(upload_to='blog/posts/', verbose_name='Gönderi Fotosu:')),
                ('title', models.CharField(max_length=250, verbose_name='Başlık:')),
                ('text', ckeditor_uploader.fields.RichTextUploadingField(verbose_name='İçerik')),
                ('slug', autoslug.fields.AutoSlugField(editable=False, populate_from='title', unique=True)),
                ('tags', models.CharField(blank=True, max_length=255)),
                ('publish', models.DateTimeField(default=django.utils.timezone.now, verbose_name='Yayınlama Tarihi')),
                ('created', models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')),
                ('updated', models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')),
                ('status', models.CharField(choices=[('DF', 'Taslak'), ('PB', 'Yayınla')], default='DF', max_length=2, verbose_name='Yayınlanma Durumu')),
                ('author', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL, verbose_name='Yazar')),
                ('category', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='post', to='post.category', verbose_name='Kategori')),
            ],
            options={
                'verbose_name': 'Gönderi',
                'verbose_name_plural': 'Gönderi',
                'db_table': 'posts',
                'ordering': ['publish'],
            },
        ),
        migrations.AddIndex(
            model_name='posts',
            index=models.Index(fields=['publish'], name='posts_publish_cbdb19_idx'),
        ),
    ]
