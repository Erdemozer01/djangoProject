from django import forms
from bioinformatic.models import Fasta, Genbank

FILE_TYPE = (
    ("fasta", "fasta"),
    ("genbank", "genbank"),
)


# creating a form
class FileTypeSelect(forms.Form):
    geeks_field = forms.ChoiceField(choices=FILE_TYPE)


class FileReadForm(forms.Form):
    file = forms.FileField()


class FastaIdForm(forms.ModelForm):
    class Meta:
        model = Fasta
        fields = ('gene',)

        labels = {
            'gene': 'Gen Bölgesi Seçiniz'
        }

        widgets = {
            'gene': forms.Select(choices=Fasta.objects.values_list('gene'))
        }


class GenbankIdForm(forms.ModelForm):
    class Meta:
        model = Genbank
        fields = ['gene']

        labels = {
            'gene': 'Gen Bölgesi Seçiniz'
        }

        widgets = {
            'gene': forms.Select(choices=Genbank.objects.all().values_list('gene', 'gene'))
        }


class XmlFileForm(forms.Form):
    file_input = forms.FileField(label="XML Dosyası Seçiniz")