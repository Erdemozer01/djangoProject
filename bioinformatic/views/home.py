from django.views import generic
from bioinformatic.models import Slide
from django.shortcuts import render


class BioinformaticHomeView(generic.ListView):
    template_name = "bioinformatic/home.html"
    model = Slide

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(BioinformaticHomeView, self).get_context_data(**kwargs)
        context['bre'] = "Biyoinformatik"
        return context
