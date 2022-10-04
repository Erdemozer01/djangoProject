import Bio
from django.shortcuts import render, redirect
from bioinformatic.forms.writing import FastaWritingForm
from bioinformatic.forms.add import AddFastaData
from django.views import generic
from Bio import SeqIO
from pathlib import Path
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from bioinformatic.forms.translation import DNAFastaFileTranslateForm

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "files", f"{f.name}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def fasta_writing(request):
    fastaform = FastaWritingForm(request.POST or None)
    if request.method == "POST":
        if fastaform.is_valid():

            id = fastaform.cleaned_data["id"]
            descriptions = fastaform.cleaned_data["description"]
            sequence = fastaform.cleaned_data["sequence"]
            sequence = Seq(sequence)
            bad_chars = [';', ':', '!', "*", "\n", '"', "\r"]

            for i in bad_chars:
                sequence = sequence.replace(i, '')

            rec1 = SeqRecord(
                sequence,
                id=id,
                description=descriptions
            )

            file = os.path.join(BASE_DIR, 'files\\file.fasta')

            SeqIO.write(rec1, file, "fasta")

            return redirect("bioinformatic:download")

        else:

            msg = "Bir hata meydana geldi"

            return render(request, 'bioinformatic/fasta/notfound.html', {
                "msg": msg
            })

    return render(request, "bioinformatic/fasta/writing.html", {
        "form": fastaform,
        "bre": "Fasta Dosyası Yazma"
    })


def append_new_line(file_name, text_to_append):
    """Append given text as a new line at the end of file"""
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        # Move read cursor to the start of file.
        file_object.seek(0)
        # If file is not empty then append '\n'
        data = file_object.read(100000)
        if len(data) > 0:
            file_object.write("\n")
        # Append text at the end of file
        file_object.write(text_to_append)


def fasta_add(request):
    form = AddFastaData(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            try:
                handle_uploaded_file(request.FILES["file"])
                input_file = form.cleaned_data['file']
                fasta_id = form.cleaned_data["fasta_id"]
                description = form.cleaned_data['description']
                sequence = form.cleaned_data["sequence"]
                file_fasta = os.path.join(BASE_DIR, "files", f"{input_file}")

                record = SeqRecord(
                    seq=Seq(sequence),
                    id=fasta_id.encode().decode(encoding="utf-8", errors="ignore"),
                    description=description
                ).format("fasta")

                read_input = open(file_fasta, "r").read()

                if "LOCUS" in read_input:
                    os.remove(file_fasta)
                    msg = "Lütfen Fasta Dosyası Seçiniz"
                    return render(request, "bioinformatic/fasta/notfound.html", {
                        "msg": msg
                    })
                if "#NEXUS" in read_input:
                    os.remove(file_fasta)
                    msg = "Lütfen Fasta Dosyası Seçiniz."
                    return render(request, "bioinformatic/fasta/notfound.html", {
                        "msg": msg
                    })
                append_new_line(file_fasta, str(record))
                os.rename(file_fasta, os.path.join(BASE_DIR, "files", "file.fasta"))
                return redirect("bioinformatic:fasta_download")

            except UnicodeDecodeError:
                os.remove(file_fasta)
                msg = "Lütfen Fasta Dosyası Seçiniz."
                return render(request, "bioinformatic/fasta/notfound.html", {
                    "msg": msg
                })

        else:

            msg = "Bir hata meydana geldi"

            return render(request, "bioinformatic/fasta/notfound.html", {
                "msg": msg
            })

    return render(request, "bioinformatic/fasta/add.html", {
        "form": form,
        "bre": "Fasta Dosyası Veri Ekleme"
    })


def fasta_file_translate(request):
    form = DNAFastaFileTranslateForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            try:
                handle_uploaded_file(request.FILES['file'])
                input_fasta_path = os.path.join(BASE_DIR, "files", f"{form.cleaned_data['file']}")
                protein_fasta_path = os.path.join(BASE_DIR, "files", "protein.txt")
                records = SeqIO.parse(input_fasta_path, format="fasta")
                table = form.cleaned_data['translate_table']

                protein_file = open(protein_fasta_path, "w")

                for record in records:
                    protein_file.write(f"{record.description}")
                    protein_file.write("\n")
                    protein_file.write("Sekans Uzunlugu: " + f"{len(record.translate().seq)}")
                    protein_file.write("\n")
                    protein_file.write(f"{record.translate(table=table).seq}")
                    protein_file.write(2 * "\n")

                os.remove(input_fasta_path)

            except Bio.BiopythonWarning:
                pass


        return redirect('bioinformatic:fasta_protein_download')

    return render(request, "bioinformatic/fasta/translate.html", {'form': form})
