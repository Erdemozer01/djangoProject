from django.db import models


class Inbox(models.Model):
    name = models.CharField(max_length=20, verbose_name='Ad Soyad:')
    email = models.EmailField(verbose_name='Email:')
    topic = models.CharField(max_length=50, verbose_name='Konu:')
    phone = models.CharField(max_length=13, verbose_name='Telefon')
    message = models.TextField(verbose_name='Mesaj')

    created = models.DateTimeField(auto_now_add=True, verbose_name='Gönderilme Tarihi')

    def __str__(self):
        return self.name

    class Meta:
        ordering = ['created']
        db_table = 'inbox'
        verbose_name = 'Gelen Mesajlar'
        verbose_name_plural = 'Gelen Mesajlar'
