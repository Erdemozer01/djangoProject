from django.db import models
from django.contrib.auth.models import User
from ckeditor_uploader.fields import RichTextUploadingField


class ReadManager(models.Manager):
    def get_queryset(self):
        return super().get_queryset() \
            .filter(status=UserMessagesModel.Status.NOT_SEEN)


class UserMessagesModel(models.Model):
    class Status(models.TextChoices):
        SEEN = 'Okundu', 'Okundu'
        NOT_SEEN = 'Okunmadı', 'Okunmadı'

    title = models.CharField(max_length=150, verbose_name="Konu:", blank=True)
    sender = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='Gönderen: ',
                                  related_name='messages_sender')
    receiver = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name='Alıcı: ',
                                    related_name='messages_receiver')
    message = RichTextUploadingField(verbose_name='Mesaj', blank=False)

    status = models.CharField(max_length=50, choices=Status.choices,
                              default=Status.NOT_SEEN, verbose_name='Görülme Durumu')

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    objects = models.Manager()
    published = ReadManager()

    def __str__(self):
        return str(self.receiver)

    class Meta:
        db_table = 'messages'
        verbose_name = 'Kullanıcı Mesajları'
        verbose_name_plural = 'Kullanıcı Mesajları'
        ordering = ['-created']