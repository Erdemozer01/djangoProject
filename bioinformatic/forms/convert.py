from django import forms
from bioinformatic.models import FileFormat

choices = FileFormat.objects.values_list('name', 'name')

molecule_type = (
    ('', '-' * 30),
    ('DNA', 'DNA'),
    ('RNA', 'RNA'),
    ('protein', 'PROTEİN'),
)


class FileConvertForm(forms.Form):
    in_file = forms.FileField(label="Dosya")
    molecule_type = forms.ChoiceField(label="Molekül Tipi", choices=molecule_type)
    in_format = forms.ChoiceField(label="Dönüştürülecek dosya formatı",
                                  choices=choices)
    out_format = forms.ChoiceField(label="Hangi Formata Dönüştürmek İstiyorsunuz ?",
                                   choices=choices)


class FileConvertModelForm(forms.ModelForm):
    class Meta:
        model = FileFormat
        fields = ['name']
