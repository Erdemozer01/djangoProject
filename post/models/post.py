from django.contrib.auth.models import User
from ckeditor_uploader.fields import RichTextUploadingField
from .category import Category
from django.db import models
from autoslug.fields import AutoSlugField
from django.utils import timezone


class PublishedManager(models.Manager):
    def get_queryset(self):
        return super().get_queryset() \
            .filter(status=Posts.Status.PUBLISHED)


class Posts(models.Model):
    class Status(models.TextChoices):
        DRAFT = 'DF', 'Taslak'
        PUBLISHED = 'PB', 'Yayınla'

    cover = models.ImageField(upload_to="blog/posts/", verbose_name="Gönderi Fotosu:")
    title = models.CharField(max_length=250, verbose_name="Başlık:", blank=False)
    text = RichTextUploadingField(verbose_name='İçerik', blank=False)
    slug = AutoSlugField(populate_from="title", unique=True)
    tags = models.CharField(max_length=255, blank=True)
    author = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='Yazar')
    category = models.ForeignKey(Category, verbose_name="Kategori", related_name='post', on_delete=models.CASCADE)

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    status = models.CharField(max_length=2, choices=Status.choices,
                              default=Status.DRAFT, verbose_name='Yayınlanma Durumu')

    objects = models.Manager()  # The default manager.
    published = PublishedManager()  # Our custom manager.

    def __str__(self):
        return self.title

    class Meta:
        db_table = 'posts'
        ordering = ['publish']
        indexes = [
            models.Index(fields=['publish']),
        ]
        verbose_name_plural = 'Gönderi'
        verbose_name = 'Gönderi'


