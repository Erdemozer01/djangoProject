import os
import stat

from django.shortcuts import *
from reportlab.lib import colors
from pathlib import Path
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from bioinformatic.forms.diagram import GenomeDiagramForm, RestrictionModelForms, RestrictionModelFormset

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "bioinformatic", "files", f"{f.name}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def genome_diagram(request):
    global file_path
    diagram_form = GenomeDiagramForm(request.POST or None, request.FILES or None)
    restric_form = RestrictionModelFormset(request.POST or None)
    if request.method == "POST":
        if diagram_form.is_valid():

            try:

                handle_uploaded_file(request.FILES['file'])
                file = diagram_form.cleaned_data['file']
                sigil = diagram_form.cleaned_data['diagram_shape']
                format = diagram_form.cleaned_data['diagram_format']
                fragments = diagram_form.cleaned_data['fragment']

                file_path = os.path.join(BASE_DIR, "bioinformatic", "files", f"{file}")

                records = SeqIO.read(file_path, "genbank")

                gd_diagram = GenomeDiagram.Diagram(records.description)
                gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
                gd_feature_set = gd_track_for_features.new_set()

                for feature in records.features:
                    if feature.type != "gene":
                        # Exclude this feature
                        continue
                    if len(gd_feature_set) % 2 == 0:
                        color = colors.tomato
                    else:
                        color = colors.lightblue
                    gd_feature_set.add_feature(feature, color=color, label=True, sigil=sigil)

                gd_diagram.draw(
                    format=format,
                    orientation="landscape",
                    pagesize="A4",
                    fragments=fragments,
                    start=0,
                    end=len(records),
                    circle_core=0.5
                )

                gd_diagram.write("plasmid_circular_nice.pdf", "PDF")
                gd_diagram.write("plasmid_linear.png", "png")

            except ValueError:
                os.remove(file_path)
                msg = "Çoklu Genbank Kaydı"
                return render(request, 'bioinformatic/fasta/notfound.html',
                              {"msg": msg, 'bre': 'Hata', 'url': reverse('bioinformatic:genome_diagram')})

    return render(request, "bioinformatic/diagram/analiz.html",
                  {'diagram_form': diagram_form, 'restric_form': restric_form})


def add_extra_enzim(request):
    return render(request, "bioinformatic/diagram/add_enzym.html", {'form': RestrictionModelFormset()})
