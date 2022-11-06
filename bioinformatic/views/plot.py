from django.shortcuts import *
from django.views import generic
from django.contrib.auth.decorators import login_required
from bioinformatic.forms.plot import PlotForm
from pathlib import Path
import os
from Bio import SeqIO
from Bio.SeqUtils import GC
import pylab
from django.core.files import File
from bioinformatic.models import GraphicModels

BASE_DIR = Path(__file__).resolve().parent.parent


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "files", f"{f}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


@login_required
def plot(request):
    form = PlotForm(request.POST or None, request.FILES or None)
    obj = GraphicModels()
    if request.method == "POST":
        if form.is_valid():
            try:
                handle_uploaded_file(request.FILES['file'])
                file_format = form.cleaned_data['file_format']
                plot_type = form.cleaned_data['plot_type']
                file_path = os.path.join(BASE_DIR, "files", f"{request.FILES['file']}")
                image_path = os.path.join(BASE_DIR, "files", "{}_{}.png".format(request.user, plot_type))
                handle = open(file_path)
                records = SeqIO.parse(handle, format=file_format)

                if GraphicModels.objects.filter(user=request.user).exists():
                    GraphicModels.objects.filter(user=request.user).delete()

                if plot_type == "histogram":
                    sizes = [len(rec) for rec in records]
                    pylab.hist(sizes, bins=20)
                    pylab.title(
                        "Sekans Uzunluk Histogram"
                    )

                    pylab.xlabel("Sekans Uzunluğu (bp)")
                    pylab.ylabel("Sayı")
                    pylab.savefig(image_path)

                elif plot_type == "gc":

                    gc_values = sorted(GC(record.seq) for record in records)

                    if gc_values == []:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Hatalı Dosya Türü", 'url': reverse('bioinformatic:plot')})

                    pylab.plot(gc_values)
                    pylab.title(
                        "%GC Plot"
                    )
                    pylab.xlabel("Gen Sayısı")
                    pylab.ylabel("%GC")
                    pylab.savefig(image_path)

                elif plot_type == "dot":
                    rec_one = next(records)
                    rec_two = next(records)
                    window = 7
                    seq_one = rec_one.seq.upper()
                    seq_two = rec_two.seq.upper()
                    data = [
                        [
                            (seq_one[i: i + window] != seq_two[j: j + window])
                            for j in range(len(seq_one) - window)
                        ]
                        for i in range(len(seq_two) - window)
                    ]

                    pylab.gray()
                    pylab.imshow(data)
                    pylab.xlabel("%s (length %i bp)" % (rec_one.id, len(rec_one)))
                    pylab.ylabel("%s (length %i bp)" % (rec_two.id, len(rec_two)))
                    pylab.title("Dot Plot")
                    pylab.savefig(image_path)

                with Path(image_path).open('rb') as image_obj:
                    obj.user = request.user
                    obj.graph_type = plot_type
                    obj.format = file_format
                    if plot_type == "histogram":
                        obj.histogram_plot = File(image_obj, name="{}_plot.png".format(plot_type))
                    elif plot_type == "gc":
                        obj.gc_plot = File(image_obj, name="{}_plot.png".format(plot_type))
                    elif plot_type == "dot":
                        obj.dot_plot = File(image_obj, name="{}_plot.png".format(plot_type))
                    obj.save()

                handle.close()
                os.remove(file_path)
                image_handle = open(image_path)
                image_handle.close()
                os.remove(image_path)

                return HttpResponseRedirect(
                    reverse('bioinformatic:plot_results',
                            args=(obj.graph_type.lower(), obj.user, obj.created.date(), obj.pk)))
            except StopIteration:
                return render(request, "bioinformatic/fasta/notfound.html",
                              {'msg': "Hatalı Dosya Türü", 'url': reverse('bioinformatic:plot')})
            except RuntimeError:
                return render(request, "bioinformatic/fasta/notfound.html",
                              {'msg': "Beklenmedik hata meydana geldi", 'url': reverse('bioinformatic:plot')})
        else:
            return redirect('bioinformatic:plot')

    return render(request, "bioinformatic/plot/input.html", {'form': form})


class PlotDetailView(generic.DetailView):
    template_name = "bioinformatic/plot/results.html"
    model = GraphicModels

    def get_context_data(self, **kwargs):
        context = super(PlotDetailView, self).get_context_data(**kwargs)
        context['bre'] = "Plot Sonuçları"
        context['object'] = GraphicModels.objects.filter(user=self.request.user).latest('created')
        return context
