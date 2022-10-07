import os.path
import pandas as pd
from Bio import SeqIO
from Bio import motifs
from pathlib import Path
from Bio.Seq import Seq
from django.shortcuts import redirect, render, reverse
from Bio.motifs.jaspar import Motif
from bioinformatic.models import FastaDNAMotifModel
from bioinformatic.forms.motif import FastaDNASeqMotifForm
from django.core.files import File
from django.contrib.auth.decorators import login_required
BASE_DIR = Path(__file__).resolve().parent.parent


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "files", f"{f.name}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)

@login_required
def fasta_motif(request):
    form = FastaDNASeqMotifForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():

            handle_uploaded_file(request.FILES['file'])
            input_file = os.path.join(BASE_DIR, "files", f"{form.cleaned_data['file']}")
            motif_file = os.path.join(BASE_DIR, "files", "motif.txt")
            motif_file_path = Path(motif_file)
            matrix_position_file = os.path.join(BASE_DIR, "files", "nmp.xlsx")
            matrix_position_file_path = Path(matrix_position_file)
            specific_scoring_matrices = os.path.join(BASE_DIR, "files", "pssm.xlsx")
            specific_scoring_matrices_path = Path(specific_scoring_matrices)
            position_weight_matrices = os.path.join(BASE_DIR, "files", "pwm.xlsx")
            position_weight_matrices_path = Path(position_weight_matrices)

            seq_len = form.cleaned_data['length']

            reading = SeqIO.parse(input_file, "fasta")

            instances = []

            for i in reading:
                instances.append(Seq(f"{i.seq[:seq_len]}"))

            motif = motifs.create(instances=instances)

            with open(motif_file, "w") as motif_obj:
                motif_obj.write("Motif : " + "\n")
                motif_obj.write(str(motif))
                motif_obj.write("\n")
                motif_obj.write("Consensus Sequence            :  ")
                motif_obj.write(str(motif.consensus))
                motif_obj.write("\n")
                motif_obj.write("Anticonsensus Sequence        :  ")
                motif_obj.write(str(motif.anticonsensus))
                motif_obj.write("\n")
                motif_obj.write("Degenerate Consensus Sequence :  ")
                motif_obj.write(str(motif.degenerate_consensus))

            pd.DataFrame(motif.counts).to_excel(matrix_position_file)
            pd.DataFrame(motif.pssm).to_excel(specific_scoring_matrices)
            pd.DataFrame(motif.pwm).to_excel(position_weight_matrices)

            if FastaDNAMotifModel.objects.all().filter(user=request.user.id).exists():
                FastaDNAMotifModel.objects.all().filter(user=request.user.id).all().delete()

            doc = FastaDNAMotifModel()

            with motif_file_path.open('r') as motif_file_obj:
                doc.motif_file = File(motif_file_obj, name=motif_file_path.name)
                doc.user = request.user
                doc.save()

            with matrix_position_file_path.open('rb') as matrix_file_obj:
                doc.mpf = File(matrix_file_obj, name=matrix_position_file_path.name)
                doc.save()

            with specific_scoring_matrices_path.open('rb') as pssm_file_obj:
                doc.pssm = File(pssm_file_obj, name=specific_scoring_matrices_path.name)
                doc.save()

            with position_weight_matrices_path.open('rb') as pwm_file_obj:
                doc.pwm = File(pwm_file_obj, name=position_weight_matrices_path.name)
                doc.save()

            results = FastaDNAMotifModel.objects.all().filter(user=request.user.id).latest('created')

            return render(request, "bioinformatic/motif/results.html", {'results': results, 'bre': "Fasta DNA Motif"})

        else:
            form = FastaDNASeqMotifForm()

    return render(request, "bioinformatic/motif/dna_seq_motif.html", {'form': form, 'bre': "Fasta DNA Motif"})
