from django import forms

FILE_TYPE = (
    ('fasta', 'FASTA'),
    ('genbank', 'GENBANK'),
)

PLOT = (
    ('gc', "%GC PLOT"),
    ('histogram', "HİSTOGRAM PLOT"),
    ('dot', "DOT PLOT"),
)

class PlotForm(forms.Form):
    file = forms.FileField(label="Dosya")
    file_type = forms.ChoiceField(label="Dosya Türü", choices=FILE_TYPE)
    plot_type = forms.ChoiceField(label="Plot Türü", choices=PLOT)
