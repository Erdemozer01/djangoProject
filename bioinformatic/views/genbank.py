from django.shortcuts import render, redirect
from bioinformatic.forms.file import FileReadForm, GenbankIdForm
from bioinformatic.forms.writing import FastaWritingForm
from bioinformatic.forms.add import AddFastaData
from Bio import SeqIO
from bioinformatic.models import Genbank
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


def genbank_read(request):
    form = FileReadForm(request.POST or None, request.FILES or None)
    if request.method == "POST":

        if form.is_valid():

            file = os.path.join(BASE_DIR, 'files\\{}'.format(form.cleaned_data['file']))
            handle_uploaded_file(request.FILES['file'])
            records = SeqIO.parse(file, "genbank")

            for record in records:
                Genbank.objects.create(gene=record.id, sekans=record.seq, description=record.description,
                                       dbxrefs=record.dbxrefs, features=record.features)

            open(file, 'r').close()

            os.remove(file)

            if Genbank.objects.exists():
                return redirect('bioinformatic:genbank_region')

            else:
                return render(request, 'bioinformatic/genbank/notfound.html', {'msg': 'Hatalı Dosya'})

        else:

            form = FileReadForm(request.POST or None, request.FILES or None)

    return render(request, 'bioinformatic/genbank/read.html', {'form': form, 'bre': 'Genbank Dosyası Okuması'})


def delete_genbank(request):
    Genbank.objects.all().delete()
    return redirect('bioinformatic:genbank_read')


def genbank_region_find(request):
    genbank = GenbankIdForm(request.POST)
    if request.method == "POST":

        if genbank.is_valid():
            sequence = Genbank.objects.get(gene=genbank.cleaned_data['gene'])

            return render(request, 'bioinformatic/fasta/sequence.html', {'seq': sequence, 'len': len(sequence.sekans)})

    return render(request, 'bioinformatic/fasta/id.html', {'form': genbank, 'bre': 'Genbank Dosyası Okuması'})


def genbank_writing(request):
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


def append_multiple_lines(file_name, lines_to_append):
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        appendEOL = False
        # Move read cursor to the start of file.
        file_object.seek(0)
        # Check if file is not empty
        data = file_object.read(100)
        if len(data) > 0:
            appendEOL = True
        # Iterate over each string in the list
        for line in lines_to_append:
            # If file is not empty then append '\n' before first line for
            # other lines always append '\n' before appending line
            if appendEOL == True:
                file_object.write("\n")
            else:
                appendEOL = True
            # Append element at the end of file
            file_object.write(line)
            file_object.close()


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

                print(rec1)

                file = os.path.join(BASE_DIR, 'files\\{}'.format(request.FILES['file']))
                new = os.path.join(BASE_DIR, 'files\\file.fasta')

                os.rename(file, new)

                append_multiple_lines(new, rec1)

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
