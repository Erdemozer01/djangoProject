from django.db import models
from django.utils import timezone


# Create your models here.

class Social(models.Model):
    name = models.CharField(max_length=20, verbose_name='Sosyal Medya', help_text='ör: Facebook')
    url = models.URLField(verbose_name='URL')

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    def __str__(self):
        return self.name

    class Meta:
        db_table = 'social'
        verbose_name = 'Sosyal Medya'
        verbose_name_plural = 'Sosyal Medya'
