import subprocess
import sys

from django.shortcuts import *
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

BASE_DIR = Path(__file__).resolve().parent.parent.parent
path = os.path.join(BASE_DIR, "bioinformatic\\files\\")


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def MuscleTreesView(request):
    global file, muscle_exe
    form = PhyloGeneticTreeForm(request.POST or None, request.FILES or None)

    if request.method == "POST":

        if form.is_valid():

            try:

                file = os.path.join(BASE_DIR, "bioinformatic\\files\\{}".format(form.cleaned_data['files']))

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
                if sys.platform.startswith('win32'):
                    muscle_exe = os.path.join(BASE_DIR, 'bioinformatic/apps/muscle3.8.425_win32.exe')
                elif sys.platform.startswith('linux'):
                    muscle_exe = os.path.join(BASE_DIR, 'bioinformatic/apps/muscle3.8.31_i86linux32')

                input_file = os.path.join(BASE_DIR, 'bioinformatic/files/{}'.format(form.cleaned_data['files']))
                output_file = os.path.join(BASE_DIR, 'bioinformatic/files/aligned.fasta')
                align_file = os.path.join(BASE_DIR, 'bioinformatic/files/align.txt')
                tree_file = os.path.join(BASE_DIR, "bioinformatic/files/tree.xml")

                muscle_cline = MuscleCommandline(muscle_exe, input=input_file, out=output_file)

                AlignIO.convert(output_file, "fasta", align_file, "clustal")

                reading_align = open(align_file, "r")

                alignment = AlignIO.read(reading_align, "clustal")

                calculator = DistanceCalculator('identity')

                constructor = DistanceTreeConstructor(calculator, method=method)

                input_file = open(input_file)

                input_file.close()

                os.remove(file)

                reading_align.close()

                os.remove(output_file)

                reading_align.close()

                os.remove(align_file)

                tree = constructor.build_tree(alignment)

                tree.rooted = True

                tree.root.color = "#596CFF"

                Phylo.write(tree, tree_file, "phyloxml")

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

                os.remove(tree_file)

                return render(request, "bioinformatic/trees/result.html",
                              {"bre": "Filogenetik Ağaç"})

            except UnicodeDecodeError:
                os.remove(file)

                return render(request, 'bioinformatic/fasta/notfound.html', {'msg': 'Hatalı Dosya Seçtiniz'})

    return render(request, "bioinformatic/trees/muscle.html",
                  {'form': form, "bre": "Muscle Filogenetik Ağaç Oluşturma"})
