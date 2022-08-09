import os.path
from pathlib import Path
from Bio import pairwise2
from Bio.Align import substitution_matrices
from django.shortcuts import render, redirect

from bioinformatic \
    .forms.alignments import GlobalForm, LocalForm

BASE_DIR = Path(__file__).resolve().parent.parent

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

            filepath = os.path.join(BASE_DIR, 'files\\global_alignment.txt')
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
            filepath = os.path.join(BASE_DIR, 'files\\local_alignment.txt')
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
