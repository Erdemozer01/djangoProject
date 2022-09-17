import datetime
import io
import os.path
import sys
from pathlib import Path
import Bio.Application
from matplotlib import pyplot as plt
from Bio import pairwise2, Phylo, SeqIO
from Bio.Align import substitution_matrices
from django.shortcuts import render, redirect, reverse
from Bio import AlignIO
import subprocess
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Align.Applications import MuscleCommandline, ClustalwCommandline, ClustalOmegaCommandline
from bioinformatic \
    .forms.alignments import GlobalForm, LocalForm, MultipleSequenceAlignmentForm, MultipleFileReading

BASE_DIR = Path(__file__).resolve().parent.parent.parent
path = os.path.join(BASE_DIR, 'bioinformatic', 'files/')


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def global_alignment(request):
    form = GlobalForm(request.POST or None)
    if request.POST:
        if form.is_valid():
            seq1 = form.cleaned_data["seq1"]
            seq2 = form.cleaned_data["seq2"]

            bad_chars = [';', ':', '!', "*", "\n", '"', "\r", " "]
            for i in bad_chars:
                seq1 = seq1.replace(i, '')
                seq2 = seq2.replace(i, '')

            alignments = pairwise2.align.globalxx(seq1, seq2)
            result = pairwise2.format_alignment(*alignments[0], full_sequences=True)

            filepath = os.path.join(BASE_DIR, 'bioinformatic\\files\\global_alignment.txt')
            file = open(filepath, "w")

            file.write(result)

            file.close()

            return redirect("bioinformatic:global_alignment_download")

        else:
            form = GlobalForm()

    return render(request, "bioinformatic/alignments/global.html", {"form": form, "bre": "Global Alignment"})


def local_alignment(request):
    form = LocalForm(request.POST or None)
    if request.POST:
        if form.is_valid():
            seq1 = form.cleaned_data["seq1"]
            seq2 = form.cleaned_data["seq2"]
            align = form.cleaned_data["algo"]

            bad_chars = [';', ':', '!', "*", "\n", '"', "\r", " "]
            for i in bad_chars:
                seq1 = seq1.replace(i, '')
                seq2 = seq2.replace(i, '')

            blosum62 = substitution_matrices.load(f"{align}")
            alignments = pairwise2.align.localds(seq1, seq2, blosum62, -10, -1)
            result = pairwise2.format_alignment(*alignments[0], full_sequences=True)
            filepath = os.path.join(BASE_DIR, 'bioinformatic\\files\\local_alignment.txt')
            file = open(filepath, "w")

            file.write(f"{align} Matrisine göre dizileme yapılmıştır\n\n")
            file.close()

            app = open(filepath, "a")
            app.write(result)

            app.close()
            return redirect("bioinformatic:local_alignment_download")

        else:
            form = LocalForm()

    return render(request, "bioinformatic/alignments/local.html", {"form": form, "bre": "Local Alignment"})


def MultipleSeqAlignment(request):
    global clustalw2_exe, muscle_exe, clustal_omega_exe, clustal_result, clustalw2
    form = MultipleSequenceAlignmentForm(request.POST or None, request.FILES or None)
    if request.method == "POST":
        if form.is_valid():

            handle_uploaded_file(request.FILES['file'])
            method = form.cleaned_data['method']
            algoritma = form.cleaned_data['algoritma']

            if method == "MUSCLE":
                try:
                    if sys.platform.startswith('win32'):
                        muscle_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'muscle3.8.425_win32.exe')
                    elif sys.platform.startswith('linux'):
                        muscle_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'muscle3.8.425_i86linux32')

                    input_file = os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                              '{}'.format(form.cleaned_data['file']))
                    output_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta')
                    align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "align.aln")
                    tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')

                    records = SeqIO.parse(input_file, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiplesequence_alignments')})

                    muscle_result = subprocess.check_output([muscle_exe, "-in", input_file, "-out", output_file])

                    AlignIO.convert(output_file, 'fasta', align_file, 'clustal')
                    alignment = AlignIO.read(align_file, "clustal")
                    calculator = DistanceCalculator('identity')

                    constructor = DistanceTreeConstructor(calculator, method=algoritma)
                    tree = constructor.build_tree(alignment)
                    Phylo.write(tree, tree_file, "phyloxml")

                    Phylo.draw(tree, do_show=False)

                    plt.xlabel('Dal uzunluğu')
                    plt.ylabel('Taksonomi')
                    if algoritma == "nj":
                        plt.title('Neighbor Joining Ağacı')
                    elif algoritma == "upgma":
                        plt.title('UPGMA Ağacı')

                    plt.suptitle(f'{method} Metodu')

                    plt.savefig(os.path.join(BASE_DIR, "media", "tree.jpg"))

                    os.remove(input_file)
                    os.remove(output_file)

                    with open(align_file, 'a') as file_obj:
                        file_obj.write('Muscle Metodu ile oluşturulmuştur.\n')
                        file_obj.write('Tarih:')
                        from django.utils import timezone
                        file_obj.write(str(timezone.now().date()))

                    return render(request, 'bioinformatic/alignments/muscle.html',
                                  {'bre': 'Muscle Metodu Sonuçları'})

                except Bio.Application.ApplicationError:
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiplesequence_alignments')})

            elif method == "clustalw2":
                try:

                    if sys.platform.startswith('win32'):
                        clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2.exe')
                    elif sys.platform.startswith('linux'):
                        clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2')

                    input_file = os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                              '{}'.format(form.cleaned_data['file']))
                    output_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta')
                    align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "align.aln")
                    dnd_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "tree.dnd")
                    tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')
                    stats = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'stats.txt')
                    pwmatrix = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'pwmatrix.txt')
                    scores_path = os.path.join(BASE_DIR, "bioinformatic", "files", "scrores.txt")


                    records = SeqIO.parse(input_file, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiplesequence_alignments')})

                    clustalw_cline = ClustalwCommandline(
                        clustalw2_exe,
                        infile=input_file,
                        outfile=output_file,
                        align=True,
                        outorder="ALIGNED",
                        convert=True,
                        output="FASTA",
                        newtree=dnd_file,
                        stats=stats
                    )

                    assert os.path.isfile(clustalw2_exe), "Clustal W executable missing"
                    stdout, stderr = clustalw_cline()

                    AlignIO.convert(output_file, 'fasta', align_file, 'clustal')
                    alignment = AlignIO.read(align_file, 'clustal')

                    calculator = DistanceCalculator('identity')
                    constructor = DistanceTreeConstructor(calculator, method=algoritma)
                    tree = constructor.build_tree(alignment)

                    Phylo.write(tree, tree_file, "phyloxml")

                    Phylo.draw(tree, do_show=False)

                    plt.xlabel('Dal uzunluğu')
                    plt.ylabel('Taksonomi')

                    if algoritma == "nj":
                        plt.title('Neighbor Joining Ağacı')
                    elif algoritma == "upgma":
                        plt.title('UPGMA Ağacı')

                    plt.suptitle(f'{method.upper()} Metodu')

                    plt.savefig(os.path.join(BASE_DIR, "media", "tree.jpg"))

                    os.remove(input_file)
                    os.remove(tree_file)
                    os.remove(dnd_file)
                    os.remove(output_file)
                    with open(align_file, 'a') as file_obj:
                        file_obj.write('Clustalw2 Metodu ile oluşturulmuştur.\n')
                        file_obj.write('Tarih:')
                        from django.utils import timezone
                        file_obj.write(str(timezone.now().date()))

                    read_stats = open(stats, 'r').readlines()[1:16]

                    os.remove(stats)

                    for i in read_stats:
                        open(stats, 'a').writelines(i)

                    open(scores_path, "w").writelines(stdout)

                    scores = open(scores_path, 'r').readlines()[44:52]

                    os.remove(scores_path)

                    for i in scores:
                        open(scores_path, 'a').writelines(i)


                    return render(request, 'bioinformatic/alignments/clustal.html',
                                  {'bre': 'Clustalw2 Metodu Sonuçları'})

                except Bio.Application.ApplicationError:
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiplesequence_alignments')})

            elif method == "omega":
                try:

                    if sys.platform.startswith('win32'):
                        clustal_omega_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps',
                                                         'clustal-omega-1.2.2-win64/clustalo.exe')
                    elif sys.platform.startswith('linux'):
                        clustal_omega_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps',
                                                         'clustalo-1.2.4-Ubuntu-32-bit')

                    input_file = os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                              '{}'.format(form.cleaned_data['file']))
                    output_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta')
                    align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "align.aln")
                    tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')

                    records = SeqIO.parse(input_file, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiplesequence_alignments')})

                    clustal_omega_cline = ClustalOmegaCommandline(clustal_omega_exe, infile=input_file,
                                                                  outfile=output_file, version=True)
                    if sys.platform.startswith('win32'):
                        assert os.path.isfile(clustal_omega_exe)
                        stdout, stderr = clustal_omega_cline()
                    elif sys.platform.startswith('linux'):
                        subprocess.Popen(str(clustal_omega_cline), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, universal_newlines=True,
                                         shell=(sys.platform != "win32"))

                    AlignIO.convert(output_file, 'fasta', align_file, 'clustal')
                    alignment = AlignIO.read(align_file, "clustal")
                    calculator = DistanceCalculator('identity')

                    constructor = DistanceTreeConstructor(calculator, method=algoritma)
                    tree = constructor.build_tree(alignment)

                    Phylo.write(tree, tree_file, "phyloxml")

                    Phylo.draw(tree, do_show=False)

                    plt.xlabel('Dal uzunluğu')
                    plt.ylabel('Taksonomi')
                    if algoritma == "nj":
                        plt.title('Neighbor Joining Ağacı')
                    elif algoritma == "upgma":
                        plt.title('UPGMA Ağacı')

                    plt.suptitle(f'{method.upper()} Metodu')

                    plt.savefig(os.path.join(BASE_DIR, "media", "tree.jpg"))

                    os.remove(input_file)
                    os.remove(tree_file)
                    with open(align_file, 'a') as file_obj:
                        file_obj.write('OMEGA Metodu ile oluşturulmuştur.\n')
                        file_obj.write('Tarih:')
                        from django.utils import timezone
                        file_obj.write(str(timezone.now().date()))

                    return render(request, 'bioinformatic/alignments/omega.html', {'bre': 'Omega Metodu Sonuçları'})

                except Bio.Application.ApplicationError:
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiplesequence_alignments')})

    return render(request, 'bioinformatic/alignments/multiple.html', {'form': form, 'bre': 'Multiple Sekans Alignment'})
