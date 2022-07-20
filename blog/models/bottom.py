from django.db import models
from django.utils import timezone


class Bottom(models.Model):
    title = models.CharField(max_length=100, verbose_name='Başlık', help_text='Anasayfa')
    text = models.TextField(verbose_name='İçerik')

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    def __str__(self):
        return self.title

    class Meta:
        ordering = ['publish']
        db_table = 'bottom_cover'
        verbose_name = 'Alt Kapak'
        verbose_name_plural = 'Alt Kapak'
