import os.path
import sys
from django.http import HttpResponseRedirect
from django.core.files import File
from pathlib import Path
import Bio.Application
from matplotlib import pyplot as plt
from Bio import pairwise2, Phylo, SeqIO
from Bio.Align import substitution_matrices
from django.shortcuts import render, redirect, reverse
from Bio import AlignIO
import subprocess
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Align.Applications import ClustalwCommandline, ClustalOmegaCommandline
from bioinformatic \
    .forms.alignments import GlobalForm, LocalForm, MultipleSequenceAlignmentForm
from bioinformatic.models import MultipleSequenceAlignment, MaximumFileSize
from django.contrib.auth.decorators import login_required
import dash_cytoscape as cyto
from dash import Dash, html, Input, Output
from django.urls import reverse_lazy

BASE_DIR = Path(__file__).resolve().parent.parent.parent
path = os.path.join(BASE_DIR, 'bioinformatic', 'files/')
msa_path = os.path.join(BASE_DIR, 'media', 'msa')


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def global_alignment(request):
    if MaximumFileSize.objects.exists():
        max_file_size = MaximumFileSize.objects.all().latest('created')
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
            from Bio import Align
            aligner = Align.PairwiseAligner()
            alignments2 = aligner.align(seq1, seq2)

            filepath = os.path.join(BASE_DIR, 'bioinformatic\\files\\global_alignment.txt')
            file = open(filepath, "w")

            file.write(result)
            file.write("\n\n")
            file.write(str(alignments2[0].substitutions))

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

            matrice = substitution_matrices.load(f"{align}")
            alignments = pairwise2.align.localds(seq1, seq2, matrice, -10, -1)
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


@login_required
def MultipleSeqAlignment(request):
    global clustalw2_exe, muscle_exe, clustal_omega_exe, clustal_result, clustalw2, clustalx_exe
    obj = MultipleSequenceAlignment()
    form = MultipleSequenceAlignmentForm(request.POST or None, request.FILES or None)

    if request.method == "POST":
        if form.is_valid():
            handle_uploaded_file(request.FILES['file'])
            method = form.cleaned_data['method']
            tree_type = form.cleaned_data['tree_type']
            molecule_type = form.cleaned_data['molecule_type']
            alignment_filetype = form.cleaned_data['alignment_filetype']

            if method == "MUSCLE":
                try:
                    user_path = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user))
                    if Path(user_path).exists():
                        pass
                    else:
                        os.makedirs(os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user)))

                    if sys.platform.startswith('win32'):
                        muscle_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'muscle3.8.425_win32.exe')
                    elif sys.platform.startswith('linux'):
                        muscle_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'muscle3.8.425_i86linux32')

                    input_file = os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                              '{}'.format(request.FILES['file']))

                    output_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligment.fasta')

                    align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', f'aligned.{alignment_filetype}')

                    if alignment_filetype == "phylip":
                        align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                                  f'aligned.phy')

                    path = Path(align_file)

                    tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')

                    if MultipleSequenceAlignment.objects.all().filter(user=request.user.id).exists():
                        MultipleSequenceAlignment.objects.all().filter(user=request.user.id).all().delete()

                    records = SeqIO.parse(input_file, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiple_sequence_alignments')})

                    muscle_result = subprocess.check_output([muscle_exe, "-in", input_file, "-out", output_file])

                    AlignIO.convert(output_file, 'fasta', align_file, f'{alignment_filetype}',
                                    molecule_type=molecule_type)
                    alignment = AlignIO.read(align_file, f'{alignment_filetype}')


                    calculator = DistanceCalculator('identity')

                    constructor = DistanceTreeConstructor(calculator, method=tree_type)
                    tree = constructor.build_tree(alignment)

                    Phylo.write(tree, tree_file, "phyloxml")

                    Phylo.draw(tree, do_show=False)

                    f = Phylo.to_networkx(tree)

                    plt.xlabel('Dal uzunluğu')
                    plt.ylabel('Taksonomi')

                    if tree_type == "nj":
                        plt.title('Neighbor Joining Ağacı')
                    elif tree_type == "upgma":
                        plt.title('UPGMA Ağacı')

                    plt.suptitle(f'{method}')

                    plt.savefig(os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                             "{}_filogenetik_ağaç.jpg".format(request.user)))

                    with path.open(mode='r') as f:
                        obj.align_file = File(f, name=path.name)
                        obj.user = request.user
                        obj.method = method
                        obj.tree_type = tree_type
                        obj.molecule_type = molecule_type
                        obj.alignment_filetype = alignment_filetype
                        obj.tree = os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                "{}_filogenetik_ağaç.jpg".format(request.user))
                        obj.save()

                    tree_path = Path(tree_file)

                    with tree_path.open(mode='rb') as tree_file_obj:
                        obj.tree_file = File(tree_file_obj, name=tree_path.name)
                        obj.save()

                    handle = open(input_file)
                    handle.close()
                    os.remove(input_file)
                    handle = open(output_file)
                    handle.close()
                    os.remove(output_file)
                    handle = open(align_file)
                    handle.close()
                    os.remove(align_file)
                    handle = open(tree_file)
                    handle.close()
                    os.remove(tree_file)

                except Bio.Application.ApplicationError:
                    os.remove(
                        os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiple_sequence_alignments')})

            elif method == "clustalw2":

                try:

                    user_path = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user))

                    if Path(user_path).exists():
                        pass
                    else:
                        os.makedirs(os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user)))

                    if sys.platform.startswith('win32'):
                        clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2.exe')
                    elif sys.platform.startswith('linux'):
                        clustalw2_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalw2')

                    input_file = os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                              '{}'.format(form.cleaned_data['file']))
                    output_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligment.fasta')

                    align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', f'aligned.{alignment_filetype}')
                    aligned_path = Path(align_file)
                    dnd_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "tree.dnd")
                    tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')
                    stats = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'stats.txt')
                    stats_path = Path(stats)
                    scores_path = os.path.join(BASE_DIR, "bioinformatic", "files", "scores.txt")
                    score_path = Path(scores_path)

                    if MultipleSequenceAlignment.objects.all().filter(user=request.user.id).exists():
                        MultipleSequenceAlignment.objects.all().filter(user=request.user.id).all().delete()

                    records = SeqIO.parse(input_file, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiple_sequence_alignments')})

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

                    AlignIO.convert(output_file, 'fasta', align_file, f'{alignment_filetype}',
                                    molecule_type=molecule_type)
                    alignment = AlignIO.read(align_file, f'{alignment_filetype}')

                    calculator = DistanceCalculator('identity')
                    constructor = DistanceTreeConstructor(calculator, method=tree_type)
                    tree = constructor.build_tree(alignment)

                    Phylo.write(tree, tree_file, "phyloxml")

                    Phylo.draw(tree, do_show=False)

                    plt.xlabel('Dal uzunluğu')
                    plt.ylabel('Taksonomi')

                    if tree_type == "nj":
                        plt.title('Neighbor Joining Ağacı')
                    elif tree_type == "upgma":
                        plt.title('UPGMA Ağacı')

                    plt.suptitle(f'{method.upper()}')

                    plt.savefig(os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                             "{}_filogenetik_ağaç.jpg".format(request.user)))

                    read_stats = open(stats, 'r').readlines()[1:16]

                    os.remove(stats)

                    for i in read_stats:
                        open(stats, 'a').writelines(i)

                    open(scores_path, "w").writelines(stdout)

                    scores = open(scores_path, 'r').readlines()[44:52]

                    os.remove(scores_path)

                    for i in scores:
                        open(scores_path, 'a').writelines(i)

                    with aligned_path.open(mode='r') as f:
                        obj.align_file = File(f, name=aligned_path.name)
                        obj.user = request.user
                        obj.method = method
                        obj.tree_type = tree_type
                        obj.molecule_type = molecule_type
                        obj.alignment_filetype = alignment_filetype
                        obj.tree = os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                "{}_filogenetik_ağaç.jpg".format(request.user))

                        obj.save()

                    with stats_path.open(mode='r') as stats_file_obj:
                        obj.stats = File(stats_file_obj, name=stats_path.name)
                        obj.save()

                    with score_path.open(mode='r') as file_obj:
                        obj.scores = File(file_obj, name=score_path.name)
                        obj.save()

                    with Path(tree_file).open(mode='r') as file_obj:
                        obj.tree_file = File(file_obj, name="tree.xml")
                        obj.save()

                    os.remove(input_file)
                    os.remove(output_file)
                    os.remove(align_file)
                    os.remove(stats)
                    os.remove(scores_path)
                    os.remove(dnd_file)
                    os.remove(tree_file)


                except Bio.Application.ApplicationError:
                    os.remove(
                        os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiple_sequence_alignments')})

            elif method == "clustalx":

                try:

                    user_path = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user))

                    if Path(user_path).exists():
                        pass
                    else:
                        os.makedirs(os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user)))

                    if sys.platform.startswith('win32'):
                        clustalx_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalx.exe')
                    elif sys.platform.startswith('linux'):
                        clustalx_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps', 'clustalx')

                    input_file = os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                              '{}'.format(form.cleaned_data['file']))
                    output_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligment.fasta')

                    align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', f'aligned.{alignment_filetype}')
                    aligned_path = Path(align_file)
                    dnd_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "tree.dnd")
                    tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')
                    stats = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'stats.txt')
                    stats_path = Path(stats)
                    scores_path = os.path.join(BASE_DIR, "bioinformatic", "files", "scores.txt")
                    score_path = Path(scores_path)

                    if MultipleSequenceAlignment.objects.all().filter(user=request.user.id).exists():
                        MultipleSequenceAlignment.objects.all().filter(user=request.user.id).all().delete()

                    records = SeqIO.parse(input_file, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiple_sequence_alignments')})

                    clustalx_cline = ClustalwCommandline(
                        clustalx_exe,
                        infile=input_file,
                        outfile=output_file,
                        align=True,
                        outorder="ALIGNED",
                        convert=True,
                        output="FASTA",
                        newtree=dnd_file,
                        stats=stats,
                    )

                    assert os.path.isfile(clustalx_exe), "Clustal W executable missing"
                    stdout, stderr = clustalx_cline()

                    AlignIO.convert(output_file, 'fasta', align_file, f'{alignment_filetype}',
                                    molecule_type=molecule_type)
                    alignment = AlignIO.read(align_file, f'{alignment_filetype}')

                    calculator = DistanceCalculator('identity')
                    constructor = DistanceTreeConstructor(calculator, method=tree_type)
                    tree = constructor.build_tree(alignment)

                    Phylo.write(tree, tree_file, "phyloxml")

                    Phylo.draw(tree, do_show=False)

                    plt.xlabel('Dal uzunluğu')
                    plt.ylabel('Taksonomi')

                    if tree_type == "nj":
                        plt.title('Neighbor Joining Ağacı')
                    elif tree_type == "upgma":
                        plt.title('UPGMA Ağacı')

                    plt.suptitle(f'{method.upper()}')

                    plt.savefig(os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                             "{}_filogenetik_ağaç.jpg".format(request.user)))

                    read_stats = open(stats, 'r').readlines()[1:16]

                    os.remove(stats)

                    for i in read_stats:
                        open(stats, 'a').writelines(i)

                    open(scores_path, "w").writelines(stdout)

                    scores = open(scores_path, 'r').readlines()[44:52]

                    os.remove(scores_path)

                    for i in scores:
                        open(scores_path, 'a').writelines(i)

                    with aligned_path.open(mode='r') as f:
                        obj.align_file = File(f, name=aligned_path.name)
                        obj.user = request.user
                        obj.method = method
                        obj.tree_type = tree_type
                        obj.molecule_type = molecule_type
                        obj.alignment_filetype = alignment_filetype
                        obj.tree = os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                "{}_filogenetik_ağaç.jpg".format(request.user))

                        obj.save()

                    with stats_path.open(mode='r') as stats_file_obj:
                        obj.stats = File(stats_file_obj, name=stats_path.name)
                        obj.save()

                    with score_path.open(mode='r') as file_obj:
                        obj.scores = File(file_obj, name=score_path.name)
                        obj.save()

                    with Path(tree_file).open(mode='r') as file_obj:
                        obj.tree_file = File(file_obj, name="tree.xml")
                        obj.save()

                    os.remove(input_file)
                    os.remove(output_file)
                    os.remove(align_file)
                    os.remove(stats)
                    os.remove(scores_path)
                    os.remove(dnd_file)
                    os.remove(tree_file)


                except Bio.Application.ApplicationError:
                    os.remove(
                        os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiple_sequence_alignments')})

            elif method == "omega":

                try:
                    user_path = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user))
                    if Path(user_path).exists():
                        pass
                    else:
                        os.makedirs(os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user)))

                    if sys.platform.startswith('win32'):
                        clustal_omega_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps',
                                                         'clustal-omega-1.2.2-win64/clustalo.exe')
                    elif sys.platform.startswith('linux'):
                        clustal_omega_exe = os.path.join(BASE_DIR, 'bioinformatic', 'apps',
                                                         'clustalo-1.2.4-Ubuntu-32-bit')

                    input_file = os.path.join(BASE_DIR, 'bioinformatic', 'files',
                                              '{}'.format(form.cleaned_data['file']))
                    output_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligment.fasta')
                    output_path = Path(output_file)
                    align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', f'aligned.{alignment_filetype}')
                    aligned_path = Path(align_file)
                    tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')

                    if MultipleSequenceAlignment.objects.all().filter(user=request.user.id).exists():
                        MultipleSequenceAlignment.objects.all().filter(user=request.user.id).all().delete()

                    records = SeqIO.parse(input_file, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

                    if len(seq_id) < 3:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                       'url': reverse('bioinformatic:multiple_sequence_alignments')})

                    clustal_omega_cline = ClustalOmegaCommandline(clustal_omega_exe,
                                                                  infile=input_file,
                                                                  outfile=output_path,
                                                                  force=True,
                                                                  verbose=True,
                                                                  auto=True,
                                                                  usekimura="yes")
                    if sys.platform.startswith('win32'):
                        assert os.path.isfile(clustal_omega_exe)
                        stdout, stderr = clustal_omega_cline()

                    elif sys.platform.startswith('linux'):
                        subprocess.Popen(str(clustal_omega_cline), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, universal_newlines=True,
                                         shell=(sys.platform != "win32"))

                    AlignIO.convert(output_file, 'fasta', align_file, f'{alignment_filetype}',
                                    molecule_type=molecule_type)
                    alignment = AlignIO.read(align_file, f'{alignment_filetype}')
                    calculator = DistanceCalculator('identity')

                    constructor = DistanceTreeConstructor(calculator, method=tree_type)
                    tree = constructor.build_tree(alignment)

                    Phylo.write(tree, tree_file, "phyloxml")

                    Phylo.draw(tree, do_show=False)

                    plt.xlabel('Dal uzunluğu')
                    plt.ylabel('Taksonomi')

                    if tree_type == "nj":
                        plt.title('Neighbor Joining Ağacı')
                    elif tree_type == "upgma":
                        plt.title('UPGMA Ağacı')

                    plt.suptitle(f'{method.upper()}')

                    plt.savefig(os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                             "{}_filogenetik_ağaç.jpg".format(request.user)))

                    with aligned_path.open(mode='r') as f:
                        obj.align_file = File(f, name=aligned_path.name)
                        obj.user = request.user
                        obj.method = method
                        obj.tree_type = tree_type
                        obj.molecule_type = molecule_type
                        obj.alignment_filetype = alignment_filetype
                        obj.tree = os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                "{}_filogenetik_ağaç.jpg".format(request.user))
                        obj.save()

                    os.remove(input_file)
                    os.remove(output_path)
                    os.remove(align_file)
                    os.remove(tree_file)

                except Bio.Application.ApplicationError:
                    os.remove(
                        os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiple_sequence_alignments')})

        return HttpResponseRedirect(
            reverse('bioinformatic:msa_results',
                    args=(obj.user, obj.method.lower(), obj.molecule_type.lower(), obj.pk)))

    return render(request, 'bioinformatic/alignments/multiple.html',
                  {'form': form, 'bre': 'Multiple Sekans Alignment'})


from Bio import Phylo
import numpy as np
import plotly.graph_objs as go
from django.views import generic

class MultipleSeqDetailView(generic.DetailView):
    template_name = "bioinformatic/alignments/msa_results.html"
    model = MultipleSequenceAlignment

def phylogenetic_tree(request):
    xml_file = os.path.join(BASE_DIR, "media", "msa", f"{request.user}", f"{request.user}_tree.xml")

    def get_circular_tree_data(tree, order='level', dist=1, start_angle=0, end_angle=360, start_leaf='first'):
        """Define  data needed to get the Plotly plot of a circular tree
        """
        # tree:  an instance of Bio.Phylo.Newick.Tree or Bio.Phylo.PhyloXML.Phylogeny
        # order: tree  traversal method to associate polar coordinates to its nodes
        # dist:  the vertical distance between two consecutive leafs in the associated rectangular tree layout
        # start_angle:  angle in degrees representing the angle of the first leaf mapped to a circle
        # end_angle: angle in degrees representing the angle of the last leaf
        # the list of leafs mapped in anticlockwise direction onto circles can be tree.get_terminals()
        # or its reversed version tree.get_terminals()[::-1].
        # start leaf: is a keyword with two possible values"
        # 'first': to map  the leafs in the list tree.get_terminals() onto a circle,
        #         in the counter-clockwise direction
        # 'last': to map  the leafs in the  list, tree.get_terminals()[::-1]

        start_angle *= np.pi / 180  # conversion to radians
        end_angle *= np.pi / 180

        def get_radius(tree):
            node_radius = tree.depths()

            #  If the tree did not record  the branch lengths  assign  the unit branch length
            #  (ex: the case of a newick tree "(A, (B, C), (D, E))")
            if not np.count_nonzero(node_radius.values()):
                node_radius = tree.depths(unit_branch_lengths=True)
            return node_radius

        def get_vertical_position(tree):
            """
            returns a dict {clade: ycoord}, where y-coord is the cartesian y-coordinate
            of a  clade root in a rectangular phylogram

            """
            n_leafs = tree.count_terminals()  # Counts the number of tree leafs.

            # Assign y-coordinates to the tree leafs
            if start_leaf == 'first':
                node_ycoord = dict((leaf, k) for k, leaf in enumerate(tree.get_terminals()))
            elif start_leaf == 'last':
                node_ycoord = dict((leaf, k) for k, leaf in enumerate(reversed(tree.get_terminals())))
            else:
                raise ValueError("start leaf can be only 'first' or 'last'")

            def assign_ycoord(clade):  # compute the y-coord for the root of this clade
                for subclade in clade:
                    if subclade not in node_ycoord:  # if the subclade root hasn't a y-coord yet
                        assign_ycoord(subclade)
                node_ycoord[clade] = 0.5 * (node_ycoord[clade.clades[0]] + node_ycoord[clade.clades[-1]])

            if tree.root.clades:
                assign_ycoord(tree.root)
            return node_ycoord

        node_radius = get_radius(tree)
        node_ycoord = get_vertical_position(tree)
        y_vals = node_ycoord.values()
        ymin, ymax = min(y_vals), max(y_vals)
        ymin -= dist  # this dist subtraction is necessary to avoid coincidence of the  first and last leaf angle

        # when the interval  [ymin, ymax] is mapped onto [0, 2pi],

        def ycoord2theta(y):
            return start_angle + (end_angle - start_angle) * (y - ymin) / float(ymax - ymin)

        def get_points_on_lines(linetype='radial', x_left=0, x_right=0, y_right=0, y_bot=0, y_top=0):
            """
            - define the points that generate a radial branch and the circular arcs, perpendicular to that branch

            - a circular arc (angular linetype) is defined by 10 points on the segment of ends
            (x_bot, y_bot), (x_top, y_top) in the rectangular layout,
             mapped by the polar transformation into 10 points that are spline interpolated
            - returns for each linetype the lists X, Y, containing the x-coords, resp y-coords of the
            line representative points
            """

            if linetype == 'radial':
                theta = ycoord2theta(y_right)
                X = [x_left * np.cos(theta), x_right * np.cos(theta), None]
                Y = [x_left * np.sin(theta), x_right * np.sin(theta), None]

            elif linetype == 'angular':
                theta_b = ycoord2theta(y_bot)
                theta_t = ycoord2theta(y_top)
                t = np.linspace(0, 1, 10)  # 10 points that span the circular arc
                theta = (1 - t) * theta_b + t * theta_t
                X = list(x_right * np.cos(theta)) + [None]
                Y = list(x_right * np.sin(theta)) + [None]

            else:
                raise ValueError("linetype can be only 'radial' or 'angular'")

            return X, Y

        def get_line_lists(clade, x_left, xlines, ylines, xarc, yarc):
            """Recursively compute the lists of points that span the tree branches"""

            # xlines, ylines  - the lists of x-coords, resp y-coords of radial edge ends
            # xarc, yarc - the lists of points generating arc segments for tree branches

            x_right = node_radius[clade]
            y_right = node_ycoord[clade]

            X, Y = get_points_on_lines(linetype='radial', x_left=x_left, x_right=x_right, y_right=y_right)

            xlines.extend(X)
            ylines.extend(Y)

            if clade.clades:

                y_top = node_ycoord[clade.clades[0]]
                y_bot = node_ycoord[clade.clades[-1]]

                X, Y = get_points_on_lines(linetype='angular', x_right=x_right, y_bot=y_bot, y_top=y_top)
                xarc.extend(X)
                yarc.extend(Y)

                # get and append the lists of points representing the  branches of the descedants
                for child in clade:
                    get_line_lists(child, x_right, xlines, ylines, xarc, yarc)

        xlines = []
        ylines = []
        xarc = []
        yarc = []
        get_line_lists(tree.root, 0, xlines, ylines, xarc, yarc)
        xnodes = []
        ynodes = []

        for clade in tree.find_clades(order='preorder'):  # it was 'level'
            theta = ycoord2theta(node_ycoord[clade])
            xnodes.append(node_radius[clade] * np.cos(theta))
            ynodes.append(node_radius[clade] * np.sin(theta))

        return xnodes, ynodes, xlines, ylines, xarc, yarc

    tree = Phylo.read(xml_file, 'phyloxml')

    traverse_order = 'preorder'

    all_clades = list(tree.find_clades(order=traverse_order))
    for k in range(len((all_clades))):
        all_clades[k].id = k

    xnodes, ynodes, xlines, ylines, xarc, yarc = get_circular_tree_data(tree, order=traverse_order, start_leaf='last')

    tooltip = []
    color = []
    for clade in tree.find_clades(order=traverse_order):
        if clade.name and clade.confidence and clade.branch_length:
            tooltip.append(f"id: {clade.id}<br>name: {clade.name}<br>branch-length: {clade.branch_length}\
                        <br>confidence: {int(clade.confidence.value)}")

            color.append[clade.confidence.value]
        elif clade.name is None and clade.branch_length is not None and clade.confidence is not None:
            color.append(clade.confidence.value)
            tooltip.append(f"id: {clade.id}<br>branch-length: {clade.branch_length}\
                        <br>confidence: {int(clade.confidence.value)}")
        elif clade.name and clade.branch_length and clade.confidence is None:
            tooltip.append(f"id: {clade.id}<br>name: {clade.name}<br>branch-length: {clade.branch_length}")
            color.append(-1)
        else:
            tooltip.append('')
            color.append(-1)

    size = [9 if c != -1 else 7 for c in color]

    pl_colorscale = [[0.0, 'rgb(10,10,150)'],  # color for leafs that haven't associated a confidence
                     [0.001, 'rgb(10,10,150)'],
                     [0.001, 'rgb(214, 47, 38)'],  # in fact the colorscale starts here
                     [0.1, 'rgb(214, 47, 38)'],
                     [0.2, 'rgb(244, 109, 67)'],
                     [0.3, 'rgb(252, 172, 96)'],
                     [0.4, 'rgb(254, 224, 139)'],
                     [0.5, 'rgb(254, 254, 189)'],
                     [0.6, 'rgb(217, 239, 139)'],
                     [0.7, 'rgb(164, 216, 105)'],
                     [0.8, 'rgb(102, 189, 99)'],
                     [0.9, 'rgb(63, 170, 89)'],
                     [1.0, 'rgb(25, 151, 79)']]

    trace_nodes = dict(type='scatter',
                       x=xnodes,
                       y=ynodes,
                       mode='markers',

                       text=tooltip,
                       hoverinfo='text',
                       opacity=1)

    trace_radial_lines = dict(type='scatter',
                              x=xlines,
                              y=ylines,
                              mode='lines',
                              line=dict(color='rgb(20,20,20)', width=1),
                              hoverinfo='none')

    trace_arcs = dict(type='scatter',
                      x=xarc,
                      y=yarc,
                      mode='lines',
                      line=dict(color='rgb(20,20,20)', width=1, shape='spline'),
                      hoverinfo='none')

    layout = dict(title='Filogenetik Ağaç',
                  font=dict(family='Balto', size=14),
                  width=700,
                  height=750,
                  autosize=True,
                  showlegend=False,
                  xaxis=dict(visible=False),
                  yaxis=dict(visible=False),
                  hovermode='closest',
                  plot_bgcolor='rgb(245,245,245)',
                  margin=dict(t=75)
                  )

    fig = go.FigureWidget(data=[trace_radial_lines, trace_arcs, trace_nodes], layout=layout)

    return render(request, "bioinformatic/alignments/tree.html", {'fig':fig})

