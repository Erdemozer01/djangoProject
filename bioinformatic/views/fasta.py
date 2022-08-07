from django.shortcuts import render, redirect
from bioinformatic.forms.file import FileReadForm, FastaIdForm
from bioinformatic.forms.writing import FastaWritingForm
from bioinformatic.forms.add import AddFastaData

from Bio import SeqIO
from bioinformatic.models import Fasta
from pathlib import Path
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def fasta_read(request):
    form = FileReadForm(request.POST or None, request.FILES or None)
    if request.method == "POST":

        if form.is_valid():

            file = os.path.join(BASE_DIR, 'files\\{}'.format(form.cleaned_data['file']))
            handle_uploaded_file(request.FILES['file'])

            try:

                with open(file, 'r') as f:
                    f.read()

            except UnicodeDecodeError:
                os.remove(file)
                return render(request, 'bioinformatic/fasta/notfound.html', {'msg': 'Zip Dosyası Seçtiniz'})

            records = SeqIO.parse(os.path.join(BASE_DIR, 'files\\{}'.format(request.FILES['file'])), "fasta")

            for record in records:
                Fasta.objects.create(gene=record.id, sekans=record.seq)

            open(os.path.join(BASE_DIR, 'files\\{}'.format(request.FILES['file'])), 'r').close()
            os.remove(os.path.join(BASE_DIR, 'files\\{}'.format(request.FILES['file'])))

            if Fasta.objects.exists():
                return redirect('bioinformatic:lokus_find')
            else:
                return render(request, 'bioinformatic/fasta/notfound.html', {'msg': 'Hatalı Dosya'})

        else:

            form = FileReadForm(request.POST or None, request.FILES or None)

    return render(request, 'bioinformatic/fasta/read.html', {'form': form, 'bre': 'Fasta Dosyası Okuması'})


def delete_fasta(request):
    Fasta.objects.all().delete()
    return redirect('bioinformatic:fasta_read')


def lokus_find(request):
    fasta = FastaIdForm(request.POST)
    if request.method == "POST":

        if fasta.is_valid():
            sequence = Fasta.objects.get(gene=fasta.cleaned_data['gene'])

            return render(request, 'bioinformatic/fasta/sequence.html', {'seq': sequence, 'len': len(sequence.sekans)})

    return render(request, 'bioinformatic/fasta/id.html', {'form': fasta, 'bre': 'Fasta Dosyası Okuması'})


def fasta_writing(request):
    fastaform = FastaWritingForm(request.POST or request.GET)
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
        data = file_object.read(100)
        if len(data) > 0:
            file_object.write("\n")
        # Append text at the end of file
        file_object.write(text_to_append)


def fasta_add(request):
    form = AddFastaData(request.POST, request.FILES)
    if request.method == "POST":
        if form.is_valid():
            try:

                handle_uploaded_file(request.FILES["file"])

                id = form.cleaned_data["id"]
                descriptions = form.cleaned_data["description"]
                sequence = form.cleaned_data["sequence"]
                sequence = Seq(sequence)

                bad_chars = [';', ':', '!', "*", "\n", '"', "\r"]

                for i in bad_chars:
                    sequence = sequence.replace(i, '')

                rec1 = SeqRecord(sequence, id=id, description=descriptions)

                path = "media/{}".format(request.FILES['file'])

                print(path)

                append_new_line(path, SeqIO.write(rec1, path, "fasta"))

                return redirect("bioinformatic:download")

            except:

                pass
        else:

            msg = "Bir hata meydana geldi"

            return render(request, "bioinformatic/fasta/notfound.html", {
                "msg": msg
            })

    return render(request, "bioinformatic/fasta/add.html", {
        "form": form,
        "bre": "Fasta Dosyası Veri Ekleme"
    })
