import os
from pathlib import Path
from django.contrib.auth.decorators import login_required
from Bio import SearchIO, SeqIO
from django.shortcuts import render, redirect, get_object_or_404
from bioinformatic.forms.file import BlastXMLForm
from bioinformatic.forms.blast import BlastResultForm
from bioinformatic.views.genbank import handle_uploaded_file
from bioinformatic.models import BlastQuery, BlastQueryResults
from Bio.Blast import NCBIWWW, NCBIXML
from django.core.files import File

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')

blast_query = None
hit_id_result = None
form = BlastResultForm(instance=blast_query)

BLAST_PROGRAM = (
    ('', '------------'),
    ('blastn', 'BLASTN'),
    ('blastp', 'BLASTP'),
    ('blastx', 'BLASTX'),
    ('tblastn', 'TBLASTN'),
    ('tblastn', 'TBLASTN'),
)

BLAST_DATABASE = (
    ('', '------------'),
    ('nr', 'BLASTN'),
    ('nt', 'BLASTP'),
)


def upload_to(instance, filename):
    return 'msa/{username}/{username}_{filename}'.format(
        username=instance.user.username, filename=filename)


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

@login_required
def xml_file(request):
    global hit_id_result
    form = BlastXMLForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():
            handle_uploaded_file(request.FILES['input_file'])
            input_file_path = os.path.join(BASE_DIR, 'files', '{}'.format(form.cleaned_data['input_file']))
            blastxml_file_path = os.path.join(BASE_DIR, 'files', 'my_blast.xml')
            hsp_file_path = os.path.join(BASE_DIR, 'files', 'hsp.txt')
            program = form.cleaned_data['program']

            record = next(SeqIO.parse(input_file_path, format="fasta"))

            result_handle = NCBIWWW.qblast(f"{program}", "nt", record.format("fasta"))

            save_file = open(blastxml_file_path, "w")

            save_file.write(result_handle.read())

            save_file.close()

            blast_qresult = SearchIO.read(blastxml_file_path, "blast-xml")

            hsp_file_path = os.path.join(BASE_DIR, 'files', 'hsp.txt')

            file_obj = open(hsp_file_path, "w")

            for hsp in blast_qresult:
                file_obj.writelines(f"{hsp}" + 3 * "\n")

            file_obj.close()

            if BlastQueryResults.objects.all().filter(user=request.user.id).exists():
                BlastQueryResults.objects.all().filter(user=request.user.id).delete()

            blast_obj = BlastQueryResults()
            blastxml_file_pathlib = Path(blastxml_file_path)
            blast_hsp_pathlib = Path(hsp_file_path)

            with blastxml_file_pathlib.open('r') as blastxml:
                blast_obj.output_file = File(blastxml, name=blastxml_file_pathlib.name)
                blast_obj.user = request.user
                blast_obj.save()

            with blast_hsp_pathlib.open('r') as hsp_file_obj:
                blast_obj.hsp_file = File(hsp_file_obj, name=blast_hsp_pathlib.name)
                blast_obj.program = program
                blast_obj.save()

            results = BlastQueryResults.objects.all().filter(user=request.user.id).latest('created')

            return render(request, 'bioinformatic/xml/result.html', {'bre': 'Blast Sonuçları', 'results':results})

        else:
            form = BlastXMLForm()

    return render(request, 'bioinformatic/xml/read.html', {'form': form, 'bre': 'BlastXML Dosyası'})


def blast_result_delete(request):
    BlastQuery.objects.all().delete()
    return redirect("bioinformatic:xml_file")
