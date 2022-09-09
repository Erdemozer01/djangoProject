from django import forms
from django.conf import settings

METHOD = (
    ('', '------------'),
    ('nj', 'Neighbor Joining'),
    ('upgma', 'UPGMA'),
)


class PhyloGeneticTreeForm(forms.Form):
    files = forms.FileField(label="Fasta Dosyası Seçiniz", widget=forms.ClearableFileInput())
    method = forms.ChoiceField(choices=METHOD, label="Metod Seçiniz")
