# Generated by Django 4.1.1 on 2022-09-21 17:53

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bioinformatic', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='multiplesequencealignment',
            name='tree',
            field=models.ImageField(upload_to='', verbose_name='Filogenetik Ağaç'),
        ),
    ]