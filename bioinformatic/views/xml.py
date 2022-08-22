import os
from pathlib import Path
from Bio import SearchIO
from django.shortcuts import render, redirect, get_object_or_404
from bioinformatic.forms.file import XmlFileForm
from bioinformatic.forms.blast import BlastResultForm
from bioinformatic.views.genbank import handle_uploaded_file
from bioinformatic.models import BlastQuery

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')

blast_query = None
hit_id_result = None
form = BlastResultForm(instance=blast_query)


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


def xml_file(request):
    global hit_id_result
    form = XmlFileForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            file = os.path.join(BASE_DIR, 'files\\{}'.format(form.cleaned_data['file_input']))
            handle_uploaded_file(request.FILES['file_input'])
            try:

                blast_qresult = SearchIO.read(file, "blast-xml")

                file_path = os.path.join(BASE_DIR, 'files\\hsp.txt')

                file_obj = open(file_path, "w")

                for hsp in blast_qresult:
                    file_obj.writelines(f"{hsp}" + 3 * "\n")

                file_obj.close()

                if BlastQuery.objects.exists():
                    BlastQuery.objects.all().delete()
                    BlastQuery.objects.create(blast=blast_qresult)

                else:
                    BlastQuery.objects.create(blast=blast_qresult)

                blast_query = get_object_or_404(BlastQuery)

                form = BlastResultForm(instance=blast_query)

                return render(request, 'bioinformatic/xml/result.html', {'form': form, 'bre': 'Sonuçlar'})

            except:
                pass

            finally:
                os.remove(file)
                BlastQuery.objects.all().delete()
        else:
            form = XmlFileForm()

    return render(request, 'bioinformatic/xml/read.html', {'form': form, 'bre': 'XML Dosyası Okuma'})


def blast_result_delete(request):
    BlastQuery.objects.all().delete()
    return redirect("bioinformatic:xml_file")
