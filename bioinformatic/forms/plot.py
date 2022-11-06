from django import forms

FILE_FORMAT = (
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
    file_format = forms.ChoiceField(label="Dosya Formatı", choices=FILE_FORMAT)
    plot_type = forms.ChoiceField(label="Plot Türü", choices=PLOT)
