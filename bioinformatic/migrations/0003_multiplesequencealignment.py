# Generated by Django 4.1.1 on 2022-09-10 21:47

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bioinformatic', '0002_alter_genbankread_protein_sequence_len'),
    ]

    operations = [
        migrations.CreateModel(
            name='MultipleSequenceAlignment',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(default='Title', max_length=100)),
                ('alignment', models.TextField()),
            ],
        ),
    ]