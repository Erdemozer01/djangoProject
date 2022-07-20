from django.db import models
from django.utils import timezone
from ckeditor_uploader.fields import RichTextUploadingField


# Create your models here.

class About(models.Model):
    name = models.CharField(max_length=100, verbose_name='Hakkımızda', default='Hakkımızda')
    about = RichTextUploadingField(verbose_name='Açıklama:')

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    def __str__(self):
        return self.name

    class Meta:
        db_table = 'about'
        verbose_name = 'Hakkımızda'
        verbose_name_plural = 'Hakkımızda'
