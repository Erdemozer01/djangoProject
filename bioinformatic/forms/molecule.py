from django import forms
from bioinformatic.models import MolecularModel


class MoleculeForm(forms.Form):
    file = forms.FileField(label="Dosya Seçiniz", help_text=".pdb veya .cif uzanlı olmalıdır")


class MultipleMoleculeForm(forms.ModelForm):
    class Meta:
        model = MolecularModel
        fields = ['in_file']
        widgets = {
            'in_file': forms.ClearableFileInput(attrs={'multiple': True}),
        }

        help_texts = {
            'in_file': ".pdb veya .cif uzantılı dosya seçiniz"
        }

        labels = {
            'in_file': "Dosya Seçiniz"
        }
