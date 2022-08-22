from django import forms
from bioinformatic.models import BlastQuery


class BlastResultForm(forms.ModelForm):
    class Meta:
        model = BlastQuery
        fields = ['blast']






