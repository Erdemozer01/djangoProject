from django.db import models
from django.contrib.auth.models import User
from ckeditor_uploader.fields import RichTextUploadingField


def avatar(instance, filename):
    # file will be uploaded to MEDIA_ROOT/beat/author/<filename>
    return 'users/avatar/{0}/{1}'.format(instance.user.profile, filename)


def cover(instance, filename):
    # file will be uploaded to MEDIA_ROOT/beat/author/<filename>
    return 'users/cover/{0}/{1}'.format(instance.user.profile, filename)


class Profile(models.Model):
    Gender = [
        ('', '------'),
        ('Erkek', 'Erkek'),
        ('Kadın', 'Kadın'),
        ('Belirtmek istemiyorum', 'Belirtmek istemiyorum'),
    ]
    user = models.OneToOneField(User, on_delete=models.CASCADE, verbose_name='Kullanıcı Adı', related_name='profile')
    first_name = models.CharField(max_length=150, verbose_name="Ad")
    last_name = models.CharField(max_length=150, verbose_name="Soyad")
    email = models.EmailField(verbose_name="Email")
    cover = models.ImageField(upload_to=cover, blank=True, verbose_name='Kapak Fotosu:')
    avatar = models.ImageField(upload_to=avatar, blank=True, verbose_name='Avatar:')
    gender = models.CharField(max_length=30, choices=Gender, default="", verbose_name='Cinsiyet', blank=True)
    location = models.CharField(max_length=100, verbose_name='Yaşadığı şehir', blank=True)
    about = models.TextField(verbose_name='Hakkımda', blank=True)
    job = models.CharField(max_length=100, verbose_name='Meslek', blank=True)
    phone = models.CharField(max_length=13, verbose_name='Telefon', blank=True, unique=True, null=True)
    birth_day = models.DateField(blank=True, null=True, verbose_name='Doğum Tarihi:')
    skils = models.CharField(max_length=200, editable=True, verbose_name='Yetenekler', blank=True)

    facebook = models.URLField(verbose_name='Facebook', blank=True)
    twitter = models.URLField(verbose_name='Twitter', blank=True)
    instagram = models.URLField(verbose_name='İnstagram', blank=True)

    created = models.DateTimeField(auto_now=True, verbose_name='Katılma Tarihi', blank=True)

    def __str__(self):
        return str(self.user)

    class Meta:
        db_table = 'profile'
        verbose_name = 'Profil'
        verbose_name_plural = 'Profil'


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
