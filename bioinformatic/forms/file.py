from django import forms
from bioinformatic.models import FastaRead, GenbankRead, FileUploadModel, BlastQueryResults

FILE_TYPE = (
    ("fasta", "fasta"),
    ("genbank", "genbank"),
)

BLAST_PROGRAM = (
    ('', '------------'),
    ('blastn', 'BLASTN'),
    ('blastp', 'BLASTP'),
    ('blastx', 'BLASTX'),
    ('tblastn', 'TBLASTN'),
    ('tblastn', 'TBLASTN'),
)

BLAST_DATABASE = (
    ('', '------------'),
    ('nr', 'nr'),
    ('nt', 'nt'),
)


# creating a form

class FileTypeSelect(forms.Form):
    geeks_field = forms.ChoiceField(choices=FILE_TYPE)


class FileReadForm(forms.Form):
    file = forms.FileField(label='Dosya Seçiniz')


class FastaIdForm(forms.Form):
    gene = forms.ModelChoiceField(queryset=FastaRead.objects.all(), label='Gen Bölgesi Seçiniz')


class GenbankIdForm(forms.Form):
    gene = forms.ModelChoiceField(queryset=GenbankRead.objects.all(), label='Gen Bölgesi Seçiniz')


class XmlIdForm(forms.Form):
    hit_id = forms.Textarea()


class BlastXMLForm(forms.Form):
    program = forms.ChoiceField(choices=BLAST_PROGRAM, label="Blast Programı Seçiniz")
    database = forms.ChoiceField(choices=BLAST_DATABASE, label="Blast Veritabanı Seçiniz")
    input_file = forms.FileField(label="Fasta Dosyası")


class FileUploadModelForm(forms.ModelForm):
    class Meta:
        model = FileUploadModel
        fields = ['file']
