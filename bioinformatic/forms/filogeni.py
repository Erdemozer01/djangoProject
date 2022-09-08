from django import forms


class PhyloGeneticTreeForm(forms.Form):
    files = forms.FileField(widget=forms.ClearableFileInput(), label="Fasta Dosyası Seçiniz")
