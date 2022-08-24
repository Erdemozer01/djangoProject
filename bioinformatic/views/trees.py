from django.shortcuts import *
from bioinformatic.forms.file import FileReadForm
from Bio import Phylo
import os
from pathlib import Path
import plotly.graph_objects as go

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def trees_draw(request):
    form = FileReadForm(request.POST or None, request.FILES or None)
    if request.method == "POST":

        if form.is_valid():

            try:

                handle_uploaded_file(request.FILES['file'])

                file = os.path.join(BASE_DIR, 'files\\{}'.format(form.cleaned_data['file']))

                if not file.endswith(".xml"):
                    msg = "Hatalı Dosya uzantısı, Lütfen xml dosyası seçiniz"
                    url = reverse("bioinformatic:filogenetik_agac")
                    return render(request, 'bioinformatic/fasta/notfound.html',
                                  {"msg": msg, 'url': url})

                trees = Phylo.read(file, "phyloxml")

                plt = Phylo.draw(trees, do_show=True)

                fig = go.Figure(data=plt)

                return render(request, "bioinformatic/trees/result.html",
                              {"bre": "Filogenetik Ağaç", "fig": fig})

            except ValueError:
                msg = "Dosyada Ağaç Bulunamadı"
                url = reverse("bioinformatic:filogenetik_agac")
                return render(request, 'bioinformatic/fasta/notfound.html',
                              {"msg": msg, 'url': url})

            finally:

                os.remove(file)

    return render(request, "bioinformatic/trees/read.html", {'form': form, "bre": "Filogenetik Ağaç Oluşturma"})
