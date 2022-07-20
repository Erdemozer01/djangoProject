from django.db import models
from django.utils import timezone


class Contact(models.Model):
    title = models.CharField(max_length=20, verbose_name='Başlık')
    email = models.EmailField(verbose_name='Email')
    phone = models.CharField(max_length=13, verbose_name='Telefon')
    address = models.TextField(verbose_name='Adres')

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    def __str__(self):
        return self.title

    class Meta:
        ordering = ['publish']
        db_table = 'contact'
        verbose_name = 'İletişim'
        verbose_name_plural = 'İletişim'
