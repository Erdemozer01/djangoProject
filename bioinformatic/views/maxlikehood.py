import os
import sys
from Bio.Application import ApplicationError
from pathlib import Path
from django.shortcuts import render, reverse
from bioinformatic.models import MultipleSequenceAlignment
from Bio import SeqIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
from matplotlib import pyplot as plt
from django.core.files import File
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from bioinformatic.forms.alignments import MaximumLikeHoodForm
from Bio.Phylo.PAML import codeml
from django.contrib.auth.decorators import login_required

BASE_DIR = Path(__file__).resolve().parent.parent.parent
path = os.path.join(BASE_DIR, 'bioinformatic', 'files/')


def handle_uploaded_file(f):
    with open(path + f.name, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


@login_required
def maxlikehood(request):
    global clustalw2_exe, cml_exe
    form = MaximumLikeHoodForm(request.POST or None, request.FILES or None)

    if request.method == "POST":
        if form.is_valid():
            handle_uploaded_file(request.FILES['file'])
            molecule_type = form.cleaned_data['molecule_type']
            tree_type = form.cleaned_data['tree_type']
            file = form.cleaned_data['file']
            palm_tools = form.cleaned_data['palm_tools']

            if palm_tools == "codeml":

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

                    align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.txt')
                    aligned_path = Path(align_file)
                    dnd_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "tree.dnd")
                    tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')
                    stats = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'stats.txt')
                    stats_path = Path(stats)
                    scores_path = os.path.join(BASE_DIR, "bioinformatic", "files", "scores.txt")
                    score_path = Path(scores_path)
                    max_like_path = os.path.join(BASE_DIR, "bioinformatic", "files", "result_{}.txt".format(palm_tools))

                    if MultipleSequenceAlignment.objects.all().filter(user=request.user.id).exists():
                        MultipleSequenceAlignment.objects.all().filter(user=request.user.id).all().delete()

                    records = SeqIO.parse(input_file, "fasta")

                    seq_id = []

                    for record in records:
                        seq_id.append(record.id)

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

                    AlignIO.convert(output_file, 'fasta', align_file, 'clustal', molecule_type=molecule_type)
                    alignment = AlignIO.read(align_file, 'clustal')

                    doc = MultipleSequenceAlignment()

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

                    plt.suptitle('Maximum Likelihood')

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

                    if palm_tools == "codeml":

                        from Bio.Phylo.PAML._paml import PamlError

                        try:
                            cml = codeml.Codeml()

                            cml.alignment = os.path.join(BASE_DIR, "bioinformatic", "files", "aligment.fasta")
                            cml.tree = os.path.join(BASE_DIR, "bioinformatic", "files", "tree.xml")
                            cml.ctl_file = os.path.join(BASE_DIR, "bioinformatic", "files", "codeml.ctl")
                            cml.out_file = os.path.join(BASE_DIR, "bioinformatic", "files",
                                                        "result_{}.txt".format(palm_tools))
                            cml.working_dir = os.path.join(BASE_DIR, "bioinformatic", "files")

                            if sys.platform.startswith('win32'):
                                cml_exe = os.path.join(BASE_DIR, "bioinformatic", "apps", "paml4.9j/bin/codeml.exe")
                                results = cml.run(verbose=True, command=cml_exe, parse=True)
                            elif sys.platform.startswith('linux'):
                                cml_exe = os.path.join(BASE_DIR, "bioinformatic", "apps", "palm_linux/codeml")
                                cml.read_ctl_file(os.path.join(BASE_DIR, "bioinformatic", "files", "codeml.ctl"))
                                stdout, stderr = cml.run(verbose=True, command=cml_exe, parse=True)


                        except PamlError:
                            pass

                    with aligned_path.open(mode='r') as f:
                        doc.align_file = File(f, name=aligned_path.name)
                        doc.user = request.user
                        doc.method = 'PALM'
                        doc.tree_type = tree_type
                        doc.molecule_type = molecule_type
                        doc.palm_tools = palm_tools
                        doc.tree = os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                "{}_filogenetik_ağaç.jpg".format(request.user))
                        doc.save()

                    with stats_path.open(mode='r') as stats_file_obj:
                        doc.stats = File(stats_file_obj, name=stats_path.name)
                        doc.save()

                    with score_path.open(mode='r') as file_obj:
                        doc.scores = File(file_obj, name=score_path.name)
                        doc.save()

                    max_like = Path(max_like_path)

                    with max_like.open(mode='r') as max_file:
                        doc.ml_file = File(max_file, name=max_like.name)
                        doc.save()

                    results = MultipleSequenceAlignment.objects.all().filter(user=request.user).latest('created')

                    return render(request, "bioinformatic/alignments/palm_results.html",
                                  {'results': results, 'bre': "Maximum Likelihood Sonuçları"})

                except ApplicationError:
                    os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))

                    return render(request, 'bioinformatic/alignments/file_notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz.', 'bre': 'Hata',
                        'url': reverse('bioinformatic:maximum_likelihood')})

                except FileNotFoundError:
                    return render(request, 'bioinformatic/fasta/notfound.html', {
                        'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                        'url': reverse('bioinformatic:multiple_sequence_alignments')})

    return render(request, "bioinformatic/alignments/palm.html", {'form': form, 'bre': 'Maximum Likelihood '})
