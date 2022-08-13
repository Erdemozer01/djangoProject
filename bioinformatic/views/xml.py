import os
from pathlib import Path
from Bio import SearchIO
from django.shortcuts import render
from bioinformatic.forms.file import XmlFileForm
from bioinformatic.views.genbank import handle_uploaded_file
import subprocess as sp
BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def xml_file(request):
    form = XmlFileForm(request.POST, request.FILES)
    if request.method == "POST":
        if form.is_valid():
            file = os.path.join(BASE_DIR, 'files\\{}'.format(form.cleaned_data['file_input']))
            handle_uploaded_file(request.FILES['file_input'])
            try:
                blast_qresult = SearchIO.read(file, "blast-xml")


                return render(request, 'bioinformatic/xml/result.html', {'result': blast_qresult})
            except:
                pass
            finally:
                os.remove(file)
        else:
            form = XmlFileForm()
    return render(request, 'bioinformatic/xml/read.html', {'form': form})
