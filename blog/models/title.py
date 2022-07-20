from django.db import models
from django.utils import timezone


# Create your models here.

class Titles(models.Model):
    name = models.CharField(max_length=20, verbose_name='Site Adı')
    type = models.CharField(max_length=20, verbose_name='Site Türü')

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    def __str__(self):
        return self.name

    class Meta:
        db_table = 'title'
        verbose_name = 'Site Başlığı'
        verbose_name_plural = 'Site Başlığı'
