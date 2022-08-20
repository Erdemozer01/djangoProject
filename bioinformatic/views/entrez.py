import os
from pathlib import Path
from django.shortcuts import render, redirect
from bioinformatic.forms.entrez import EntrezForm
from Bio import Entrez

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')


def entrez(request):
    form = EntrezForm(request.POST)

    if request.method == "POST":
        if form.is_valid():

            email = form.cleaned_data['email']
            database = form.cleaned_data['database']
            accession = form.cleaned_data['accession']
            rettype = form.cleaned_data['rettype']

            Entrez.email = f"{email}"

            handle = Entrez.efetch(db=f"{database}", id=f"{accession}", rettype=f"{rettype}", retmode="text")

            read = handle.read()

            file_path = os.path.join(BASE_DIR, 'files\\entrez.txt')

            file_obj = open(file_path, "w")

            file_obj.write(read)

            return redirect("bioinformatic:entrez_file_download")

        else:
            form = EntrezForm(request.POST)
    return render(request, "bioinformatic/entrez/efetch.html", {'form': form, 'bre': "Entrez Dosyası"})
