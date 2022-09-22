import datetime
import io
import uuid
from django.contrib import messages
import os.path
import sys
from django.contrib.auth.decorators import login_required
from django.views import generic
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
from Bio.Align.Applications import MuscleCommandline, ClustalwCommandline, ClustalOmegaCommandline
from bioinformatic \
    .forms.alignments import GlobalForm, LocalForm, MultipleSequenceAlignmentForm
from bioinformatic.models import MultipleSequenceAlignment

BASE_DIR = Path(__file__).resolve().parent.parent.parent
path = os.path.join(BASE_DIR, 'bioinformatic', 'files/')
msa_path = os.path.join(BASE_DIR, 'media', 'msa')


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
    if request.user.is_authenticated:
        form = MultipleSequenceAlignmentForm(request.POST or None, request.FILES or None)
        if request.method == "POST":
            if form.is_valid():
                handle_uploaded_file(request.FILES['file'])
                method = form.cleaned_data['method']
                algoritma = form.cleaned_data['algoritma']

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

                        align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.aln')

                        path = Path(align_file)

                        tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')

                        if MultipleSequenceAlignment.objects.all().filter(user=request.user, method=method):
                            MultipleSequenceAlignment.objects.all().filter(user=request.user, method=method).delete()

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
                        doc = MultipleSequenceAlignment()

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

                        plt.savefig(os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                 "{}_filogenetik_ağaç.jpg".format(request.user)))

                        with open(align_file, 'a') as file_obj:
                            file_obj.write('Muscle Metodu kullanilmistir\n')
                            file_obj.write('Tarih: ')
                            for created in MultipleSequenceAlignment.objects.all().filter(user=request.user,
                                                                                          method=method):
                                print(created.created)

                        with path.open(mode='r') as f:
                            doc.align_file = File(f, name=path.name)
                            doc.user = request.user
                            doc.method = method
                            doc.algoritma = algoritma
                            doc.tree = os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                    "{}_filogenetik_ağaç.jpg".format(request.user))
                            doc.save()

                        os.remove(input_file)
                        os.remove(output_file)
                        os.remove(align_file)
                        os.remove(tree_file)
                        results = MultipleSequenceAlignment.objects.all().filter(user=request.user).latest('created')

                        return render(request, 'bioinformatic/alignments/muscle.html',
                                      {'bre': 'Muscle Metodu Sonuçları', 'results': results})

                    except Bio.Application.ApplicationError:
                        os.remove(
                            os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                        os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                        return render(request, 'bioinformatic/fasta/notfound.html', {
                            'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                            'url': reverse('bioinformatic:multiplesequence_alignments')})

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
                        align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "aligned.aln")
                        aligned_path = Path(align_file)
                        dnd_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "tree.dnd")
                        tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')
                        stats = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'stats.txt')
                        stats_path = Path(stats)
                        pwmatrix = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'pwmatrix.txt')
                        scores_path = os.path.join(BASE_DIR, "bioinformatic", "files", "scores.txt")
                        score_path = Path(scores_path)

                        if MultipleSequenceAlignment.objects.all().filter(user=request.user, method=method):
                            MultipleSequenceAlignment.objects.all().filter(user=request.user, method=method).delete()

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

                        doc = MultipleSequenceAlignment()

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

                        plt.savefig(os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                 "{}_filogenetik_ağaç.jpg".format(request.user)))

                        with open(align_file, 'a') as file_obj:
                            file_obj.write('Clustalw2 Metodu kullanilmistir.\n')
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

                        with aligned_path.open(mode='r') as f:
                            doc.align_file = File(f, name=aligned_path.name)
                            doc.user = request.user
                            doc.method = method
                            doc.algoritma = algoritma
                            doc.tree = os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                    "{}_filogenetik_ağaç.jpg".format(request.user))
                            doc.save()

                        with stats_path.open(mode='r') as stats_file_obj:
                            doc.stats = File(stats_file_obj, name=stats_path.name)
                            doc.save()

                        with score_path.open(mode='r') as file_obj:
                            doc.scores = File(file_obj, name=score_path.name)
                            doc.save()

                        os.remove(input_file)
                        os.remove(tree_file)
                        os.remove(dnd_file)
                        os.remove(output_file)
                        os.remove(align_file)
                        os.remove(scores_path)
                        os.remove(stats)

                        results = MultipleSequenceAlignment.objects.all().filter(user=request.user).latest('created')
                        return render(request, "bioinformatic/alignments/clustal.html",
                                      {'results': results, 'bre': "ClustalW2 Metodu Sonuçları"})

                    except Bio.Application.ApplicationError:
                        os.remove(
                            os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                        os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                        return render(request, 'bioinformatic/fasta/notfound.html', {
                            'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                            'url': reverse('bioinformatic:multiplesequence_alignments')})

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
                        align_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', "aligned.aln")
                        aligned_path = Path(align_file)
                        tree_file = os.path.join(BASE_DIR, 'bioinformatic', 'files', 'tree.xml')

                        if MultipleSequenceAlignment.objects.all().filter(user=request.user, method=method, algoritma=algoritma):
                            MultipleSequenceAlignment.objects.all().filter(user=request.user, method=method, algoritma=algoritma).delete()

                        records = SeqIO.parse(input_file, "fasta")

                        doc = MultipleSequenceAlignment()

                        seq_id = []

                        for record in records:
                            seq_id.append(record.id)

                        if len(seq_id) < 3:
                            return render(request, "bioinformatic/fasta/notfound.html",
                                          {'msg': "Ağaç oluşturmak için en az 3 canlı türü olmalıdır.",
                                           'url': reverse('bioinformatic:multiplesequence_alignments')})

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

                        AlignIO.convert(output_file, 'fasta', aligned_path, 'clustal')
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

                        plt.savefig(os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                 "{}_filogenetik_ağaç.jpg".format(request.user)))

                        with open(align_file, 'a') as file_obj:
                            file_obj.write('OMEGA Metodu kullanilmistir.\n')
                            file_obj.write('Tarih:')
                            from django.utils import timezone
                            file_obj.write(str(timezone.now().date()))

                        with aligned_path.open(mode='r') as f:
                            doc.align_file = File(f, name=aligned_path.name)
                            doc.user = request.user
                            doc.method = method
                            doc.algoritma = algoritma
                            doc.tree = os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                                    "{}_filogenetik_ağaç.jpg".format(request.user))
                            doc.save()
                            f.close()

                            os.remove(input_file)
                            os.remove(tree_file)
                            os.remove(output_file)
                            os.remove(aligned_path)

                            results = MultipleSequenceAlignment.objects.all().filter(user=request.user, method=method, algoritma=algoritma).latest(
                                'created')

                        return render(request, 'bioinformatic/alignments/omega.html',
                                      {'bre': 'Omega Metodu Sonuçları', 'results': results})

                    except Bio.Application.ApplicationError:
                        os.remove(
                            os.path.join(BASE_DIR, 'bioinformatic', 'files', '{}'.format(form.cleaned_data['file'])))
                        os.remove(os.path.join(BASE_DIR, 'bioinformatic', 'files', 'aligned.fasta'))

                        return render(request, 'bioinformatic/fasta/notfound.html', {
                            'msg': 'Hatalı Dosya Seçtiniz. Lütfen fasta dosyası seçiniz.',
                            'url': reverse('bioinformatic:multiplesequence_alignments')})

        return render(request, 'bioinformatic/alignments/multiple.html',
                      {'form': form, 'bre': 'Multiple Sekans Alignment'})
    else:
        messages.error(request, 'Multiple Sekans Alignment için giriş yapınız')
        return redirect('login')


class MultipleSequenceAlignmentListView(generic.ListView):
    model = MultipleSequenceAlignment
    template_name = "bioinformatic/alignments/clustal.html"

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(MultipleSequenceAlignmentListView, self).get_context_data(**kwargs)
        context['laborant'] = MultipleSequenceAlignment.objects.all().filter(user=self.request.user).latest('created')

        context['bre'] = "Multiple Sequence Alignment Sonuçları"
        return context
