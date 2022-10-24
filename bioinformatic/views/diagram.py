from django.shortcuts import *
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from bioinformatic.forms.diagram import GenomeDiagramForm, RestrictionModelForms, RestrictionModelFormset


def genome_diagram(request):
    diagram_form = GenomeDiagramForm(request.POST or None, request.FILES or None)
    restric_form = RestrictionModelFormset(request.POST or None)
    return render(request, "bioinformatic/diagram/analiz.html",
                  {'diagram_form': diagram_form, 'restric_form': restric_form})


def add_extra_enzim(request):
    return render(request, "bioinformatic/diagram/add_enzym.html", {'form':RestrictionModelFormset()})