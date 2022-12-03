import plotly
from django.contrib import messages
from django.shortcuts import *
from django.views import generic
from django.contrib.auth.decorators import login_required
from bioinformatic.forms.plot import PlotForm, PlotSelectForm
from pathlib import Path
import os
from Bio import SeqIO
from Bio.SeqUtils import GC
import pylab
from django.core.files import File
from bioinformatic.models import GraphicModels
from django_plotly_dash import DjangoDash
from dash import html, dcc, ctx
import plotly.express as px
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent


def handle_uploaded_file(f):
    with open(os.path.join(BASE_DIR, "files", f"{f}"), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def plot_select(request):
    if request.user.is_anonymous:
        from django.conf import settings
        messages.error(request, "Lütfen Giriş Yapınız")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
    form = PlotSelectForm(request.POST or None)
    object = GraphicModels()
    obj = GraphicModels.objects.filter(user=request.user)

    if request.method == "POST":
        if form.is_valid():
            if obj.exists():
                obj.delete()
            graph_type = form.cleaned_data['plot_type']
            object.user = request.user
            object.graph_type = graph_type
            object.save()
            obj = GraphicModels.objects.filter(user=request.user).latest('created')
            return HttpResponseRedirect(
                reverse('bioinformatic:plot',
                        args=(obj.pk, obj.user, obj.graph_type)))
    return render(request, "bioinformatic/plot/plot.html", {'form': form, 'bre': 'Plot Türü Seçiniz'})


@login_required
def plot(request, pk, user, graph_type):
    form = PlotForm(request.POST or None, request.FILES or None)
    obj = GraphicModels.objects.filter(user=request.user).latest('created')
    if request.method == "POST":
        if form.is_valid():
            try:
                handle_uploaded_file(request.FILES['file'])
                file_format = form.cleaned_data['file_format']
                file_path = os.path.join(BASE_DIR, "files", f"{request.FILES['file']}")
                handle = open(file_path)
                records = SeqIO.parse(handle, format=file_format)
                obj.format = file_format
                obj.save()

                if obj.graph_type == "histogram":
                    seq_len = [len(rec) for rec in records]

                    df = pd.DataFrame({
                        'Sekans Uzunluğu': seq_len,
                    })

                    app = DjangoDash("histogram")
                    app.layout = html.Div([
                        dcc.Graph(
                            figure=px.histogram(df, x="Sekans Uzunluğu", title="Sekans Uzunluk Dağılımı"),
                        )
                    ])

                elif obj.graph_type == "gc":
                    app = DjangoDash("gc_plot")

                    app2 = DjangoDash('name')

                    gc_values = sorted(GC(record.seq) for record in records)

                    name = []

                    gc = []

                    for rec in SeqIO.parse(file_path, format=file_format):
                        name.append(rec.name)
                        gc.append(GC(rec.seq))

                    data = pd.DataFrame({
                        'Organizma': name,
                        '%GC': gc
                    })

                    if gc_values == []:
                        return render(request, "bioinformatic/fasta/notfound.html",
                                      {'msg': "Hatalı Dosya Türü", 'url': reverse('bioinformatic:plot')})

                    df = pd.DataFrame({
                        '%GC': gc_values,
                    })

                    app.layout = html.Div([
                        dcc.Graph(
                            figure=px.line(df, title="%GC Plot", y="%GC"),
                        )
                    ])

                    app2.layout = html.Div([
                        dcc.Graph(
                            figure=px.scatter(data, title="%GC Plot", y="%GC", color="Organizma"),
                        )
                    ])

                elif obj.graph_type == "dot":
                    app = DjangoDash("dot")

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
                    import plotly.graph_objects as go

                    app.layout = html.Div([
                        dcc.Graph(
                            figure=px.imshow(data),
                            responsive=True
                        )
                    ])

                handle.close()
                os.remove(file_path)

                return HttpResponseRedirect(
                    reverse('bioinformatic:plot_results',
                            args=(obj.graph_type, obj.user, obj.created.date(), obj.pk)))
            except StopIteration:
                return render(request, "bioinformatic/fasta/notfound.html",
                              {'msg': "Hatalı Dosya Türü", 'url': reverse('bioinformatic:plot')})
            except RuntimeError:
                return render(request, "bioinformatic/fasta/notfound.html",
                              {'msg': "Beklenmedik hata meydana geldi", 'url': reverse('bioinformatic:plot')})
        else:
            return redirect('bioinformatic:plot')

    return render(request, "bioinformatic/plot/input.html", {'form': form, 'obj': obj, 'bre': f"{obj.graph_type.upper()} Plot Oluşturma"})


class PlotDetailView(generic.DetailView):
    template_name = "bioinformatic/plot/results.html"
    model = GraphicModels

    def get_context_data(self, **kwargs):
        context = super(PlotDetailView, self).get_context_data(**kwargs)
        obj = GraphicModels.objects.filter(user=self.request.user).latest('created')
        context['bre'] = f"{obj.graph_type.title()} Plot Sonuçları"
        context['object'] = GraphicModels.objects.filter(user=self.request.user).latest('created')
        return context
