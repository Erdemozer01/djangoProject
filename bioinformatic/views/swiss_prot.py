import os
from pathlib import Path
from bioinformatic.models import SwissProtModel
import Bio.SwissProt
from django.views import generic
import tempfile
import shutil

from bioinformatic.forms.file import FileReadForm
from bioinformatic.forms.url import UrlForm
from django.shortcuts import *
from Bio import SwissProt, ExPASy
from urllib.request import urlopen

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def handle_uploaded_file(file):
    with open(path + file, 'wb+') as destination:
        for chunk in file.chunks():
            destination.write(chunk)
            destination.close()


def swiss_prot_file(request):
    form = FileReadForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            try:

                file = form.cleaned_data["file"]

                handle_uploaded_file(file=file)

                file_path = os.path.join(path, file)

                if os.path.getsize(file_path) > 300000:
                    handle = open(file_path)

                    records = SwissProt.parse(handle)

                    if SwissProtModel.objects.exists():
                        SwissProtModel.objects.all().delete()
                        for record in records:
                            SwissProtModel.objects.create(accessions=record.accessions[0],
                                                          taxonomy_id=record.taxonomy_id[0],
                                                          created=record.created[0], organism=record.organism,
                                                          sequence=record.sequence,
                                                          sequence_length=record.sequence_length,
                                                          annotation_update=record.annotation_update[0])

                    else:
                        for record in records:
                            SwissProtModel.objects.create(accessions=record.accessions[0],
                                                          taxonomy_id=record.taxonomy_id[0],
                                                          created=record.created[0], organism=record.organism,
                                                          sequence=record.sequence,
                                                          sequence_length=record.sequence_length,
                                                          annotation_update=record.annotation_update[0])

                    return render(request, "bioinformatic/swiss/result.html",
                                  {"file": records, "bre": "Dosya İçerği"})

                else:

                    handle = open(file_path)

                    records = SwissProt.read(handle)

                    authors = [ref.authors for ref in records.references]

                    title = [ref.title for ref in records.references]

                    handle.close()

                    os.remove(file_path)

                    return render(request, "bioinformatic/swiss/result.html",
                                  {"file": records, "authors": authors, "title": title, "bre": "Dosya İçerği"})

            except Bio.SwissProt.SwissProtParserError:
                return render(request, 'bioinformatic/fasta/notfound.html',
                              {"msg": "Hatalı Dosya Türü, Lütfen Swiss Prot Dosyası giriniz",
                               "url": reverse("bioinformatic:swiss_prot_file"),
                               "bre": "Hata"})

            finally:
                file_path = os.path.join(BASE_DIR, "files\\{}".format(form.cleaned_data["file"]))
                open(file_path).close()
                os.remove(file_path)

        else:
            form = FileReadForm()
            return reverse("bioinformatic:swiss_prot_file")

    return render(request, "bioinformatic/swiss/file.html", {"form": form, "bre": "Swiss-Prot Dosya Okuması"})


class SwissProtListView(generic.ListView):
    pass


def swiss_prot_url(request):
    form = UrlForm(request.POST or None)

    if request.method == "POST":
        if form.is_valid():

            url = form.cleaned_data["url"]

            handle = urlopen(url)

            swiss_file = open(path + "swiss_prot.txt", "w")

            for i in handle:
                swiss_file.writelines(str(i).replace("b'", ""))

            return redirect("bioinformatic:swiss_prot_download")

    return render(request, "bioinformatic/swiss/url.html", {"form": form, "bre": "Swiss-Prot Url Dosya Okuması"})
