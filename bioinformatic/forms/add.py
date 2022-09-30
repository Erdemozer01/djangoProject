from django import forms


class AddFastaData(forms.Form):
    file = forms.FileField(label='Dosya Seçiniz')

    id = forms.CharField(
        label="İD",
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'id numarası'
            }
        )
    )

    name = forms.CharField(
        label="NAME",
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': "NAME"
            }
        )
    )

    description = forms.CharField(
        label="Description",
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Tanım'
            }
        )
    )

    dbxrefs = forms.CharField(
        label="Dbxrefs",
        required=False,
        help_text="Boş Olabilir",
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Dbxrefs'
            }
        )
    )

    annotations = forms.CharField(
        label="Annotations",
        required=False,
        help_text="Boş Olabilir",
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'Annotations'
            }
        )
    )

    sequence = forms.CharField(
        label="Sequence",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sekans'
            }
        )
    )