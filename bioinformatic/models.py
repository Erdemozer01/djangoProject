from django.db import models


# Create your models here.


class LabSlideModel(models.Model):
    image = models.ImageField(upload_to="laboratory/slide/")
    title = models.CharField(max_length=100)
    text = models.CharField(max_length=100)
    content = models.TextField()

    def __str__(self):
        return self.title

    class Meta:
        db_table = "bioinformatic_home"
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
    organism = models.CharField(verbose_name="Organizma", max_length=1000)
    description = models.CharField(max_length=1000, verbose_name="Tanım")
    protein_id = models.CharField(max_length=255, verbose_name="protein_id", blank=True, null=True)
    taxonomy = models.CharField(max_length=1000, verbose_name="Taksonomi", null=True, blank=True)
    dna_sequence = models.TextField()
    dna_sequence_len = models.BigIntegerField(verbose_name="DNA Uzunluğu", null=True, blank=True)
    protein_sequence = models.TextField(blank=True, null=True)
    protein_sequence_len = models.BigIntegerField(verbose_name="Protein Uzunluğu", null=True, blank=True)

    def __str__(self):
        return self.organism

    class Meta:
        ordering = ["id"]
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


class SwissProtModel(models.Model):
    accessions = models.CharField(max_length=1000, verbose_name="Erişim Numarası")
    organism = models.CharField(max_length=1000, verbose_name="Organizma")
    sequence = models.TextField(verbose_name="Sekans")
    sequence_length = models.CharField(max_length=1000, verbose_name="Sekans Uzunluğu")

    def __str__(self):
        return self.accessions

    class Meta:
        db_table = "swiss-prot"
        verbose_name = "swiss-prot"
        verbose_name_plural = "swiss-prot"


class BigFileUploadModel(models.Model):
    big_file = models.FileField(verbose_name="Dosya Seçme", upload_to="laboratory/bigfile/")
    created = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return self.big_file

    class Meta:
        db_table = "bigfile"
        verbose_name = "bigfile"
        verbose_name_plural = "bigfile"


class FileUploadModel(models.Model):
    file = models.FileField(verbose_name="Dosya Seçme", upload_to="laboratory/bigfile/")
    created = models.DateTimeField(auto_now_add=True, null=True)

    def __str__(self):
        return self.file

    class Meta:
        db_table = "file"
        verbose_name = "file"
        verbose_name_plural = "files"


class MultipleSequenceAlignment(models.Model):
    title = models.CharField(max_length=100, default="Title")
    alignment = models.TextField()

    def __str__(self):
        return self.title

    class Meta:
        db_table = "MSA"
        verbose_name = "MultipleSequenceAlignment"
        verbose_name_plural = "MultipleSequenceAlignment"
