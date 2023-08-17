from django.db import models
from django.contrib.auth.models import User


class ReadManager(models.Manager):
    def get_queryset(self):
        return super().get_queryset() \
            .filter(status=AuthorMessagesModel.Status.NOT_SEEN)


class AuthorMessagesModel(models.Model):
    class Status(models.TextChoices):
        SEEN = 'Okundu', 'Okundu'
        NOT_SEEN = 'Okunmadı', 'Okunmadı'

    name = models.CharField(verbose_name="Ad-Soyad", max_length=100)

    email = models.EmailField(verbose_name="Email")

    receiver = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='Alıcı: ',
                                 related_name='author_messages_receiver')

    message = models.TextField(verbose_name='Mesaj', blank=False)

    status = models.CharField(max_length=50, choices=Status.choices,
                              default=Status.NOT_SEEN, verbose_name='Görülme Durumu')

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    objects = models.Manager()
    published = ReadManager()

    def __str__(self):
        return str(self.receiver)

    class Meta:
        db_table = 'author_messages'
        verbose_name = 'Yazar Mesajları'
        verbose_name_plural = 'Yazar Mesajları'
        ordering = ['-created']
