# Generated by Django 4.1.5 on 2023-08-15 14:48

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('accounts', '0009_alter_usermessagesmodel_status'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usermessagesmodel',
            name='receiver',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='messages_receiver', to=settings.AUTH_USER_MODEL, verbose_name='Alıcı: '),
        ),
        migrations.AlterField(
            model_name='usermessagesmodel',
            name='sender',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='messages_sender', to=settings.AUTH_USER_MODEL, verbose_name='Gönderen: '),
        ),
    ]