from django.db import models


# Create your models here.


class Slide(models.Model):
    image = models.ImageField(upload_to="labratory/")
    title = models.CharField(max_length=100)
    text = models.CharField(max_length=100)
    content = models.TextField()

    def __str__(self):
        return self.title

    class Meta:
        db_table = "slide"
        verbose_name = "Slide"
        verbose_name_plural = "Slide"


class Fasta(models.Model):
    name = models.CharField(max_length=100, default='Fasta')
    gene = models.CharField(max_length=255)
    sekans = models.TextField(default='')

    def __str__(self):
        return self.gene

    class Meta:
        db_table = "fasta"
        verbose_name = "Fasta"
        verbose_name_plural = "Fasta"


class Genbank(models.Model):
    name = models.CharField(max_length=100, default='Genbank')
    gene = models.CharField(max_length=255)
    sekans = models.TextField()
    description = models.CharField(max_length=1000)
    features = models.TextField(default='')
    dbxrefs = models.CharField(max_length=1000)

    def __str__(self):
        return self.gene

    class Meta:
        db_table = "genbank"
        verbose_name = "Genbank"
        verbose_name_plural = "Genbank"
