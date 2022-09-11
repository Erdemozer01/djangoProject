from django import forms
from bioinformatic.models import MultipleSequenceAlignment
MATRIS = (
    ('', '------------'),
    ('BENNER22', 'BENNER22'),
    ('BENNER6', 'BENNER6'),
    ('BENNER74', 'BENNER74'),
    ('BLOSUM45', 'BLOSUM45'),
    ('BLOSUM50', 'BLOSUM50'),
    ('BLOSUM62', 'BLOSUM62'),
    ('BLOSUM80', 'BLOSUM80'),
    ('BLOSUM90', 'BLOSUM90'),
    ('DAYHOFF', 'DAYHOFF'),
    ('FENG', 'FENG'),
    ('HOXD70', 'HOXD70'),
    ('JOHNSON', 'JOHNSON'),
    ('JONES', 'JONES'),
    ('LEVIN', 'LEVIN'),
    ('MCLACHLAN', 'MCLACHLAN'),
    ('MDM78', 'MDM78'),
    ('PAM250', 'PAM250'),
    ('PAM30', 'PAM30'),
    ('PAM70', 'PAM70'),
    ('RAO', 'RAO'),
    ('RISLER', 'RISLER'),
    ('SCHNEIDER', 'SCHNEIDER'),
    ('STR', 'STR'),
    ('TRANS', 'TRANS'),
)

METHOD = (
    ('', '------------'),
    ('clustalw2', 'clustalw2'.upper()),
    ('MUSCLE', 'MUSCLE'),
)

class LocalForm(forms.Form):
    algo = forms.ChoiceField(
        label='Dizileme Matrisi Seçiniz',
        choices=MATRIS,
        widget=forms.Select(
            attrs={
                'class': 'form-control',

            }
        )
    )

    seq1 = forms.CharField(
        label="Sekans1 Giriniz",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sequence1'
            }
        )
    )

    seq2 = forms.CharField(
        label="Sekans2 Giriniz",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sequence2'

            }
        )
    )


class GlobalForm(forms.Form):
    seq1 = forms.CharField(
        label="Sekans1 Giriniz",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sequence1'
            }
        )
    )

    seq2 = forms.CharField(
        label="Sekans2 Giriniz",
        widget=forms.Textarea(
            attrs={
                'class': 'form-control',
                'placeholder': 'Sequence2'

            }
        )
    )


class MultipleSequenceAlignmentForm(forms.Form):
    method = forms.ChoiceField(choices=METHOD, label="Multiple Sekans Alignment Metodu Seçiniz")
    file = forms.FileField(label="Fasta Dosyası Giriniz")


class MultipleFileReading(forms.ModelForm):
    class Meta:
        model = MultipleSequenceAlignment
        fields = ['alignment']