# Generated by Django 4.1.5 on 2023-08-15 00:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('accounts', '0007_alter_usermessagesmodel_title'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='usermessagesmodel',
            options={'ordering': ['-created'], 'verbose_name': 'Kullanıcı Mesajları', 'verbose_name_plural': 'Kullanıcı Mesajları'},
        ),
        migrations.AddField(
            model_name='usermessagesmodel',
            name='status',
            field=models.CharField(choices=[('görüldü', 'Görüldü'), ('görülmedi', 'Görülmedi')], default='görülmedi', max_length=50, verbose_name='Okunma Durumu'),
        ),
    ]
