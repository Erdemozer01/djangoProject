import os.path
import sys
from pathlib import Path
from matplotlib import pyplot as plt
from Bio import pairwise2, Phylo
from Bio.Align import substitution_matrices
from django.shortcuts import render, redirect
from Bio import AlignIO
import subprocess
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Align.Applications import MuscleCommandline, ClustalwCommandline, ClustalOmegaCommandline
from bioinformatic \
    .forms.alignments import GlobalForm, LocalForm, MultipleSequenceAlignmentForm, MultipleFileReading

BASE_DIR = Path(__file__).resolve().parent.parent.parent
path = os.path.join(BASE_DIR, "bioinformatic\\files\\")


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def global_alignment(request):
    form = GlobalForm(request.POST or None)
    if request.POST:
        if form.is_valid():
            seq1 = form.cleaned_data["seq1"]
            seq2 = form.cleaned_data["seq2"]

            bad_chars = [';', ':', '!', "*", "\n", '"', "\r", " "]
            for i in bad_chars:
                seq1 = seq1.replace(i, '')
                seq2 = seq2.replace(i, '')

            alignments = pairwise2.align.globalxx(seq1, seq2)
            result = pairwise2.format_alignment(*alignments[0], full_sequences=True)

            filepath = os.path.join(BASE_DIR, 'bioinformatic\\files\\global_alignment.txt')
            file = open(filepath, "w")

            file.write(result)

            file.close()

            return redirect("bioinformatic:global_alignment_download")

        else:
            form = GlobalForm()

    return render(request, "bioinformatic/alignments/global.html", {"form": form, "bre": "Global Alignment"})


def local_alignment(request):
    form = LocalForm(request.POST or None)
    if request.POST:
        if form.is_valid():
            seq1 = form.cleaned_data["seq1"]
            seq2 = form.cleaned_data["seq2"]
            align = form.cleaned_data["algo"]

            bad_chars = [';', ':', '!', "*", "\n", '"', "\r", " "]
            for i in bad_chars:
                seq1 = seq1.replace(i, '')
                seq2 = seq2.replace(i, '')

            blosum62 = substitution_matrices.load(f"{align}")
            alignments = pairwise2.align.localds(seq1, seq2, blosum62, -10, -1)
            result = pairwise2.format_alignment(*alignments[0], full_sequences=True)
            filepath = os.path.join(BASE_DIR, 'bioinformatic\\files\\local_alignment.txt')
            file = open(filepath, "w")

            file.write(f"{align} Matrisine göre dizileme yapılmıştır\n\n")
            file.close()

            app = open(filepath, "a")
            app.write(result)

            app.close()
            return redirect("bioinformatic:local_alignment_download")

        else:
            form = LocalForm()

    return render(request, "bioinformatic/alignments/local.html", {"form": form, "bre": "Local Alignment"})


def MultipleSeqAlignment(request):
    global muscle_exe
    form = MultipleSequenceAlignmentForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            handle_uploaded_file(request.FILES['file'])
            method = form.cleaned_data['method']
            if method == "MUSCLE":
                return redirect('bioinformatic:filogenetik_agac_fasta')

            elif method == "clustalw2":
                if sys.platform.startswith('win32'):
                    muscle_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2.exe')
                elif sys.platform.startswith('linux'):
                    muscle_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2')

                input_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file']))
                output_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta')
                dnd_file = os.path.join(BASE_DIR, "bioinformatic", "files", "turtles.fasta.dnd")

                open(output_file, 'w')

                clustalw_cline = ClustalwCommandline(os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2.exe'), infile=input_file, outfile=output_file, pim=True)
                assert os.path.isfile(os.path.join(BASE_DIR, "bioinformatic", "apps", "clustalw2"))
                stdout, stderr = clustalw_cline()


                align_file = os.path.join(BASE_DIR, 'bioinformatic\\files\\align.aln')
                AlignIO.convert(output_file, 'fasta', align_file, 'clustal')
                tree = Phylo.read(dnd_file, "newick")
                Phylo.draw(tree, branch_labels=lambda c: c.branch_length, do_show=False)

                plt.savefig(os.path.join(BASE_DIR, "media", "tree.jpg"))

                os.remove(input_file)
                os.remove(dnd_file)
                plt.close()

                return render(request, 'bioinformatic/alignments/multiple_result.html')
    return render(request, 'bioinformatic/alignments/multiple.html', {'form': form, 'bre': 'Multiple Sekans Alignment'})
