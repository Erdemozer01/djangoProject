from django import forms
from bioinformatic.models import BlastQuery, BlastHSP


class BlastResultForm(forms.ModelForm):
    class Meta:
        model = BlastQuery
        fields = ['blast']


class BlastHspForm(forms.ModelForm):
    class Meta:
        model = BlastHSP
        fields = ['hit_id']



