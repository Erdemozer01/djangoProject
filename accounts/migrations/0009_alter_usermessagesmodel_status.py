# Generated by Django 4.1.5 on 2023-08-15 02:39

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('accounts', '0008_alter_usermessagesmodel_options_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='usermessagesmodel',
            name='status',
            field=models.CharField(choices=[('Okundu', 'Okundu'), ('Okunmadı', 'Okunmadı')], default='Okunmadı', max_length=50, verbose_name='Görülme Durumu'),
        ),
    ]
