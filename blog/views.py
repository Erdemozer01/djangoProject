from django.contrib.auth.mixins import LoginRequiredMixin
from django.views.generic import *
from pip._internal.locations import user_site
from post.models import Posts
from django.contrib.auth.models import User
from django.shortcuts import render
from .models import *
from django.contrib.messages.views import SuccessMessageMixin
from .forms import ContactForm
from django.urls import reverse_lazy
from django.contrib import messages
from .models.profile import Profile


class ProfileDetailView(DetailView, LoginRequiredMixin):
    template_name = 'dashboard/profil.html'
    model = Profile

    def get_queryset(self):
        return Profile.objects.filter(user=self.request.user)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(ProfileDetailView, self).get_context_data(**kwargs)
        context['posts'] = Posts.objects.filter(author_id=self.request.user.id)
        return context


class HomeView(ListView):
    template_name = 'blog/home.html'
    model = Posts
    queryset = Posts.objects.all().order_by('-id')[:5]

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(HomeView, self).get_context_data(**kwargs)
        context['user'] = User.objects.all().filter(username=user_site)
        return context


def about(request):
    about = About.objects.all()
    return render(request, 'blog/pages/about.html', {'about': about})


def terms(request):
    terms = Terms.objects.all()
    return render(request, 'blog/pages/terms.html', {'terms': terms})


class ContactView(CreateView, SuccessMessageMixin, ListView):
    form_class = ContactForm
    template_name = 'blog/pages/contact.html'
    model = Contact
    context_object_name = 'contact'
    success_url = reverse_lazy('blog:contact')
    success_message = 'Mesajınız iletildi'

    def form_valid(self, form):
        form.save()
        messages.success(self.request, self.success_message)
        return super(ContactView, self).form_valid(form)
