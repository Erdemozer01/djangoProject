from django import forms
from django.forms.models import inlineformset_factory

DIAGRAM_SHAPE = (
    ('', '------------'),
    ('ARROW', 'OK'),
    ('BOX', 'KUTU'),
)

DIAGRAM_FORMAT = (
    ('', '------------'),
    ('linear', 'DÜZ'),
    ('circular', 'DAİRESEL'),
)

from bioinformatic.models import RestrictionModel, RestrictionUserModel


class GenomeDiagramForm(forms.Form):
    file = forms.FileField(label="Genbank Dosyası", help_text="Çoklu kayıt girmeyiniz")
    diagram_shape = forms.ChoiceField(label="Diagram Şekli", choices=DIAGRAM_SHAPE)
    diagram_format = forms.ChoiceField(label="Diagram Formatı", choices=DIAGRAM_FORMAT)
    fragment = forms.IntegerField(label="Fragment Sayısı", initial=1, min_value=1)


class RestrictionModelForms(forms.ModelForm):
    class Meta:
        model = RestrictionModel
        exclude = ("laborant",)


RestrictionModelFormset = inlineformset_factory(
    RestrictionUserModel, RestrictionModel, RestrictionModelForms, min_num=1, extra=0, can_delete=False
)
