from django import forms


class FastaDNASeqMotifForm(forms.Form):
    file = forms.FileField(label="DNA Fasta Dosyası Giriniz")
    length = forms.IntegerField(label="Motif Sekans Uzunluğu")
