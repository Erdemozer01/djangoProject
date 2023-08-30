from django.views.generic import *
from pip._internal.locations import user_site
from post.models import Posts
from django.contrib.auth.models import User
from django.shortcuts import *
from .models import *
from django.contrib.messages.views import SuccessMessageMixin
from .forms import ContactForm
from django.urls import reverse_lazy
from django.contrib import messages
from .models.profile import Profile
from .forms import AuthorMessagesForm


def blog_profile_detail(request, pk, user):
    profile = Profile.objects.get(user=pk)
    author_posts = Posts.objects.filter(author=pk)
    form = AuthorMessagesForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            title = form.cleaned_data['title']
            message = form.cleaned_data['message']
            email = form.cleaned_data['email']
            sender = request.user
            receiver = User.objects.get(username=user)
            if str(request.user.username) == str(user):
                messages.error(request, 'Kendinize mesaj gönderdiniz.')
                return redirect(request.META['HTTP_REFERER'])

            AuthorMessagesModel.objects.create(title=title, sender=sender, receiver=receiver, message=message, email=email)
            messages.success(request, 'Mesajnız gönderilmiştir...')
            return redirect(request.META['HTTP_REFERER'])

        else:
            form = AuthorMessagesForm()

    return render(request, 'blog/pages/profil.html', {'profile': profile, 'author_posts': author_posts, 'form': form})


class ProfileDetailView(DetailView, CreateView):
    template_name = 'blog/pages/profil.html'
    model = Profile
    form_class = AuthorMessagesForm

    def form_valid(self, form):
        title = form.cleaned_data['title']
        sender = self.request.user
        receiver = User.objects.get(username=self.request.user)
        message = form.cleaned_data['message']
        AuthorMessagesModel.objects.create(title=title, sender=sender, receiver=receiver, message=message)
        messages.success(self.request, 'Mesajnız gönderilmiştir...')
        return super().form_valid(form)

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super(self).get(request, *args, **kwargs)

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

