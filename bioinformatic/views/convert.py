from pathlib import Path
import os
from django.shortcuts import *
from bioinformatic.forms.convert import FileConvertForm, AddFileFormatModelForm
from Bio import SeqIO

BASE_DIR = Path(__file__).resolve().parent.parent.parent


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "bioinformatic", "files", f"{f}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def file_convert(request):
    form = FileConvertForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            handle_uploaded_file(request.FILES['in_file'])
            in_format = form.cleaned_data['in_format']
            in_file_path = os.path.join(BASE_DIR, "bioinformatic", "files", f"{request.FILES['in_file']}")
            out_file_path = os.path.join(BASE_DIR, "bioinformatic", "files", "converted_{}".format(in_format))
            in_format = form.cleaned_data['in_format']
            out_format = form.cleaned_data['out_format']
            molecule_type = form.cleaned_data['molecule_type']
            SeqIO.convert(
                in_file=in_file_path,
                in_format=in_format,
                out_file=out_file_path,
                out_format=out_format,
                molecule_type=molecule_type
            )

    return render(request, "bioinformatic/convert/input.html", {'form': form})


from django.contrib.auth.decorators import permission_required, login_required
from bioinformatic.models import FileFormat
from django.contrib.auth import settings

@login_required
@permission_required('bioinformatic.add_file_format', login_url=settings.LOGIN_URL, raise_exception=True)
def add_file_format(request):
    form = AddFileFormatModelForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            file_format = form.cleaned_data['name']
            if FileFormat.objects.filter(name=file_format).exists():
                pass
            else:
                form.save()
            return redirect('bioinformatic:file_convert')
        else:
            form = AddFileFormatModelForm()

    return render(request, "bioinformatic/convert/add_format.html", {'form': form})
