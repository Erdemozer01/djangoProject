from django.db import models
from django.utils import timezone
from ckeditor_uploader.fields import RichTextUploadingField


class Terms(models.Model):
    name = models.CharField(max_length=100, verbose_name='Koşul ve Şartlar')
    terms = RichTextUploadingField(verbose_name='Koşul ve Şartlar')

    publish = models.DateTimeField(default=timezone.now, verbose_name='Gönderilme Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['created']
        verbose_name = "Koşul ve Şartlar"
        verbose_name_plural = "Koşul ve Şartlar"
