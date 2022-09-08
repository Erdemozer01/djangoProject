from django.shortcuts import *
from bioinformatic.forms.file import FileReadForm
from Bio import Phylo, SeqIO
import os
from pathlib import Path
import matplotlib.pyplot as plt
from django.conf import settings
from Bio import AlignIO
from django.views import generic
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from bioinformatic.forms.filogeni import PhyloGeneticTreeForm
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')

muscle_exe = os.path.join(settings.MUSCLE_DIR)


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def trees_draw(request):
    form = FileReadForm(request.POST or None, request.FILES or None)
    file = None
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

                tree = Phylo.read(file, "phyloxml")

                Phylo.draw(tree, do_show=False)

                plt.xlabel('Dal uzunluğu')
                plt.ylabel('Taksonomi')

                img_path = os.path.join(settings.MEDIA_ROOT, "tree.jpg")

                plt.savefig(img_path)

                plt.show()

                return render(request, "bioinformatic/trees/result.html",
                              {"bre": "Filogenetik Ağaç"})

            except ValueError:
                msg = "Dosyada Ağaç Bulunamadı"
                url = reverse("bioinformatic:filogenetik_agac")
                return render(request, 'bioinformatic/fasta/notfound.html',
                              {"msg": msg, 'url': url})

            finally:

                os.remove(file)

    return render(request, "bioinformatic/trees/xml.html",
                  {'form': form, "bre": "XML Dosyasından Filogenetik Ağaç Oluşturma"})


def FastaCreateTreesView(request):
    global file, reading_align
    form = PhyloGeneticTreeForm(request.POST or None, request.FILES or None)

    if request.method == "POST":

        if form.is_valid():

            try:

                file = os.path.join(BASE_DIR, 'files\\{}'.format(form.cleaned_data['files']))

                handle_uploaded_file(form.cleaned_data['files'])

                if not ".fasta" in file:
                    msg = "Hatalı Dosya uzantısı, Lütfen .fasta uzantılı dosyası seçiniz"
                    url = reverse("bioinformatic:filogenetik_agac_fasta")
                    return render(request, 'bioinformatic/fasta/notfound.html',
                                  {"msg": msg, 'url': url})

                records = SeqIO.parse(file, "fasta")

                seq_id = []

                sequence = []

                for record in records:
                    seq_id.append(record.id)

                if len(seq_id) < 3:
                    return render(request, "bioinformatic/fasta/notfound.html", {'msg': "Ağaç oluşturmak"
                                                                                        " için en az 3 canlı türü olmalıdır.",
                                                                                 'url': reverse(
                                                                                     'bioinformatic:filogenetik_agac_fasta')})
                muscle_cline = MuscleCommandline(muscle_exe, input=file, out=path + "aligned.fasta")

                muscle_cline()

                AlignIO.convert(path + "aligned.fasta", "fasta", path + "aligned.aln", "clustal")

                os.remove(path + "aligned.fasta")

                reading_align = open(path + "aligned.aln", "r")

                alignment = AlignIO.read(reading_align, "clustal")

                calculator = DistanceCalculator('identity')

                constructor = DistanceTreeConstructor(calculator)

                tree = constructor.build_tree(alignment)

                tree.rooted = True

                Phylo.write(tree, path + "tree.xml", "phyloxml")

                Phylo.draw(tree, do_show=False)

                img_path = os.path.join(settings.MEDIA_ROOT, "tree.jpg")

                plt.savefig(img_path)

                return render(request, "bioinformatic/trees/result.html",
                              {"bre": "Filogenetik Ağaç"})

            except UnicodeDecodeError:
                os.remove(file)

                return render(request, 'bioinformatic/fasta/notfound.html', {'msg': 'Hatalı Dosya Seçtiniz'})



    return render(request, "bioinformatic/trees/fasta.html",
                  {'form': form, "bre": "Fasta Dosyasından Filogenetik Ağaç Oluşturma"})
