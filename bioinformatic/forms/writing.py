from django import forms


class FastaWritingForm(forms.Form):
    id = forms.CharField(
        label="id",
        widget=forms.TextInput(
            attrs={
                'class': 'form-control',
                'placeholder': 'id numarası'
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

    sequence = forms.CharField(
        label="Sequence",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sekans'
            }
        )
    )
