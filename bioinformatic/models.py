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


class FastaRead(models.Model):
    gene = models.CharField(max_length=1000)
    sequence = models.TextField()

    def __str__(self):
        return self.gene

    class Meta:
        db_table = "fasta_read"
        verbose_name = "Fasta Okuması"
        verbose_name_plural = "Fasta Okumaları"


class GenbankRead(models.Model):
    gene = models.CharField(max_length=255)
    sequence = models.TextField()
    description = models.CharField(max_length=1000)
    features = models.TextField()
    dbxrefs = models.CharField(max_length=1000)

    def __str__(self):
        return self.gene

    class Meta:
        db_table = "genbank_read"
        verbose_name = "Genbank Okuması"
        verbose_name_plural = "Genbank Okumaları"


class BlastQuery(models.Model):
    blast = models.TextField()

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    class Meta:
        ordering = ['-created']


class PubMedArticle(models.Model):
    email = models.EmailField()

    article_id = models.CharField(max_length=100)

    title = models.CharField(max_length=1000)

    link = models.URLField(default="https://pubmed.ncbi.nlm.nih.gov/")

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.link

    class Meta:
        verbose_name = "PubMed Makale"
        verbose_name_plural = "PubMed Makale"


class MedlineArticle(models.Model):
    email = models.EmailField()

    article_id = models.CharField(max_length=1000)

    doi = models.CharField(max_length=1000, default="")

    title = models.CharField(max_length=1000)

    term = models.CharField(max_length=1000)

    created = models.DateTimeField(auto_now_add=True, verbose_name='Oluşturulma Tarihi')

    def __str__(self):
        return self.title





