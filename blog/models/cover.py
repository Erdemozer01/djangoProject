from django.db import models
from django.utils import timezone
from django.utils.translation import gettext_lazy as _


class Cover(models.Model):
    THEME = [
        ('', '------'),
        ('dark', _('Siyah')),
        ('primary', _('Pembe')),
        ('danger', _('Kırmızı')),
        ('info', _('Mavi')),
        ('success', _('Yaşil')),
    ]
    image = models.ImageField(upload_to='home/cover/', verbose_name='Anasayfa Kapak', blank=True)
    title = models.CharField(max_length=100, verbose_name='Başlık', help_text='Anasayfa')
    text = models.TextField(verbose_name='İçerik')
    theme = models.CharField(verbose_name=_('Tema'), choices=THEME, max_length=100, default='', blank=True,
                             help_text=_('Kapak Fotoğrafı seçtiyseniz seçin'))

    publish = models.DateTimeField(default=timezone.now, verbose_name='Yayınlama Tarihi')
    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')
    updated = models.DateTimeField(auto_now=True, verbose_name='Güncellenme Tarihi')

    def __str__(self):
        return self.title

    class Meta:
        ordering = ['publish']
        db_table = 'cover'
        verbose_name = 'Üst Kapak'
        verbose_name_plural = 'Üst Kapak'
