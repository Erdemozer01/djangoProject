from django import forms
from bioinformatic.models import FastaRead, GenbankRead, BlastHSP

FILE_TYPE = (
    ("fasta", "fasta"),
    ("genbank", "genbank"),
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

class XmlFileForm(forms.Form):
    file_input = forms.FileField(label="XML Dosyası Seçiniz")
