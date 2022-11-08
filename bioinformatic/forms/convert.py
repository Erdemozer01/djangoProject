from django import forms
from bioinformatic.models import FileFormat

file_formats = FileFormat.objects.values_list('name', 'name')

molecule_type = (
    ('', '-'*30),
    ('DNA', 'DNA'),
    ('RNA', 'RNA'),
    ('protein', 'PROTEİN'),
)


class FileConvertForm(forms.Form):
    in_file = forms.FileField(label="Dosya")
    molecule_type = forms.ChoiceField(label="Molekül Tipi", choices=molecule_type)
    in_format = forms.ChoiceField(label="Dönüştürülecek dosya formatı", choices=file_formats)
    out_format = forms.ChoiceField(label="Hangi Formata Dönüştürmek İstiyorsunuz ?", choices=file_formats)


class AddFileFormatModelForm(forms.ModelForm):
    class Meta:
        model = FileFormat
        fields = ['name']