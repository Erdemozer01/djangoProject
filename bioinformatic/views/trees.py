import subprocess

from django.shortcuts import *
from bioinformatic.forms.file import FileReadForm
from Bio import Phylo, SeqIO
import os
from pathlib import Path
import matplotlib.pyplot as plt
from django.conf import settings
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from bioinformatic.forms.filogeni import PhyloGeneticTreeForm
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files')

muscle_exe = os.path.join(settings.MUSCLE_DIR)


def handle_uploaded_file(f):
    with open(path + "/" + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def trees_draw(request):
    form = FileReadForm(request.POST or None, request.FILES or None)
    file = None
    if request.method == "POST":

        if form.is_valid():

            try:

                handle_uploaded_file(request.FILES['file'])

                file = os.path.join(BASE_DIR, 'files'.format(form.cleaned_data['file']))

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

                file = os.path.join(path, '{}'.format(form.cleaned_data['files']))

                handle_uploaded_file(form.cleaned_data['files'])

                method = form.cleaned_data['method']

                if not ".fasta" in file:
                    msg = "Hatalı Dosya uzantısı, Lütfen .fasta uzantılı dosyası seçiniz"
                    url = reverse("bioinformatic:filogenetik_agac_fasta")
                    return render(request, 'bioinformatic/fasta/notfound.html',
                                  {"msg": msg, 'url': url})

                records = SeqIO.parse(file, "fasta")

                seq_id = []

                for record in records:
                    seq_id.append(record.id)

                if len(seq_id) < 3:
                    return render(request, "bioinformatic/fasta/notfound.html", {'msg': "Ağaç oluşturmak"
                                                                                        " için en az 3 canlı türü olmalıdır.",
                                                                                 'url': reverse(
                                                                                     'bioinformatic:filogenetik_agac_fasta')})

                aligned_fasta = os.path.join(path, 'aligned.fasta')

                align_file = os.path.join(path, 'aligned.aln')

                open(aligned_fasta, 'w')

                muscle_cline = MuscleCommandline(muscle_exe, input=file, out=aligned_fasta)

                import subprocess
                muscle_result = subprocess.check_output([muscle_exe, "-in", file, "-out", aligned_fasta])

                AlignIO.convert(aligned_fasta, "fasta", align_file, "clustal")

                reading_align = open(align_file, "r")

                alignment = AlignIO.read(reading_align, "clustal")

                calculator = DistanceCalculator('identity')

                constructor = DistanceTreeConstructor(calculator, method=method)

                input_file = open(path + "/" + "{}".format(form.cleaned_data['files']))

                input_file.close()

                os.remove(file)

                reading_align.close()

                os.remove(aligned_fasta)

                reading_align.close()

                os.remove(align_file)

                tree = constructor.build_tree(alignment)

                Phylo.write(tree, path + "tree.xml", "phyloxml")

                Phylo.draw(tree, do_show=False)

                plt.xlabel('Dal uzunluğu')

                plt.ylabel('Taksonomi')

                title = form.cleaned_data['method']

                if title == "nj":
                    plt.title("Neighbor Joining Ağacı")
                else:
                    plt.title(f"{title.upper()} Ağacı")

                img_path = os.path.join(settings.MEDIA_ROOT, "tree.jpg")

                plt.savefig(img_path)

                os.remove(path + "tree.xml")

                return render(request, "bioinformatic/trees/fasta_result.html",
                              {"bre": "Filogenetik Ağaç", 'muscle_result': muscle_result})

            except UnicodeDecodeError:
                os.remove(file)

                return render(request, 'bioinformatic/fasta/notfound.html', {'msg': 'Hatalı Dosya Seçtiniz'})

    return render(request, "bioinformatic/trees/fasta.html",
                  {'form': form, "bre": "Fasta Dosyasından Filogenetik Ağaç Oluşturma"})
