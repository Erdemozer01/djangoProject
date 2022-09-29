from blog.models import Inbox, Terms, Contact, Bottom, About, Cover, Social, Titles
from django.views import generic
from django.urls import reverse_lazy
from django.shortcuts import render, reverse, redirect
from django.views.generic import ListView
from post.models import Posts
from django.contrib.auth.models import User
from django.contrib import messages
from .forms import UserRegistrationForm, UserEditForm, UserCreationForm, ProfileEditForm, UserProfileEditForm
from django.contrib.auth.mixins import LoginRequiredMixin
from .models import Profile
from django.contrib.auth.decorators import login_required
from django.contrib.auth.views import PasswordChangeView
from django.contrib.auth.forms import PasswordChangeForm
from django.views.generic import DetailView, CreateView


def blog_dashboard(request):
    return render(request, 'dashboard/blog.html', context={
        'bre': "Blog Modelleri",
        'bottom': Bottom.objects.all(),
        'about': About.objects.all(),
        'terms': Terms.objects.all(),
        'contact': Contact.objects.all(),
        'title': Titles.objects.all(),
        'social': Social.objects.all(),
        'cover': Cover.objects.all()
    })

class AddContactView(CreateView):
    template_name = 'accounts/contact.html'
    model = Contact
    fields = '__all__'
    success_url = reverse_lazy('blog_dashboard')

    def get_context_data(self, **kwargs):
        context = super(AddContactView, self).get_context_data(**kwargs)
        context['bre'] = 'Site Başlığı'
        return context


def delete_contact(request):
    Contact.objects.all().delete()
    if Contact.DoesNotExist:
        return redirect('blog_dashboard')
    return redirect('blog_dashboard')


class AddCoverView(CreateView):
    template_name = 'accounts/topcover.html'
    model = Cover
    fields = '__all__'
    success_url = reverse_lazy('blog_dashboard')

    def get_context_data(self, **kwargs):
        context = super(AddCoverView, self).get_context_data(**kwargs)
        context['bre'] = 'Site Başlığı'
        return context


def delete_cover(request):
    Cover.objects.all().delete()
    if Cover.DoesNotExist:
        return redirect('blog_dashboard')
    return redirect('blog_dashboard')


class AddSocialView(CreateView):
    template_name = 'accounts/socialmedia.html'
    model = Social
    fields = '__all__'
    success_url = reverse_lazy('blog_dashboard')

    def get_context_data(self, **kwargs):
        context = super(AddSocialView, self).get_context_data(**kwargs)
        context['bre'] = 'Site Başlığı'
        return context


def delete_social(request):
    Social.objects.all().delete()
    if Social.DoesNotExist:
        return redirect('blog_dashboard')
    return redirect('blog_dashboard')


class AddTitlesView(CreateView):
    template_name = 'accounts/title.html'
    model = Titles
    fields = '__all__'
    success_url = reverse_lazy('blog_dashboard')

    def get_context_data(self, **kwargs):
        context = super(AddTitlesView, self).get_context_data(**kwargs)
        context['bre'] = 'Site Başlığı'
        return context


def delete_titles(request):
    Titles.objects.all().delete()
    if Titles.DoesNotExist:
        return redirect('blog_dashboard')
    return redirect('blog_dashboard')


class AddTermsView(CreateView):
    template_name = 'accounts/terms.html'
    model = Terms
    fields = '__all__'
    success_url = reverse_lazy('blog_dashboard')


def delete_terms(request):
    Terms.objects.all().delete()
    if Terms.DoesNotExist:
        return redirect('blog_dashboard')
    return redirect('blog_dashboard')


class AddAboutView(CreateView):
    template_name = 'accounts/about.html'
    model = About
    fields = '__all__'
    success_url = reverse_lazy('blog_dashboard')


def delete_about(request):
    About.objects.all().delete()
    if About.DoesNotExist:
        return redirect('blog_dashboard')
    return redirect('blog_dashboard')


class AddBottomView(CreateView):
    template_name = 'accounts/bottom.html'
    model = Bottom
    fields = '__all__'
    success_url = reverse_lazy('blog_dashboard')


def delete_bottom(request):
    Bottom.objects.all().delete()
    if Bottom.DoesNotExist:
        return redirect('blog_dashboard')
    return redirect('blog_dashboard')


class PostsDashBoardView(ListView, LoginRequiredMixin):
    template_name = 'accounts/posts.html'
    model = Posts
    queryset = Posts.objects.all().order_by("-id")
    login_url = reverse_lazy('login')

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(PostsDashBoardView, self).get_context_data(**kwargs)
        context['bre'] = "Gönderi Ekle"
        return context


class UsersView(ListView, LoginRequiredMixin):
    template_name = 'index.html'
    model = Posts
    queryset = Posts.objects.all().order_by("-id")

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(UsersView, self).get_context_data(**kwargs)
        context['users'] = User.objects.all()
        context['messages'] = Inbox.objects.all()
        return context


class ProfileDetailView(DetailView, LoginRequiredMixin):
    template_name = 'dashboard/profil.html'
    model = Profile
    queryset = Profile.objects.all()

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(ProfileDetailView, self).get_context_data(**kwargs)
        context['posts'] = Posts.objects.filter(author_id=self.request.path.split('/')[4])

        return context


class UserRegister(generic.CreateView):
    template_name = "accounts/register.html"
    form_class = UserRegistrationForm
    success_url = reverse_lazy("register")

    def form_valid(self, form):
        messages.success(self.request, 'Başarılı Bir Şekilde Kayıt Oldunuz.')
        form.save()
        return super(UserRegister, self).form_valid(form)


class MessageDetail(generic.DetailView):
    template_name = "pages/mesage.html"
    model = Inbox
    queryset = Inbox.objects.all()
    context_object_name = "msg"


class MessageDeleteView(generic.DeleteView):
    template_name = "pages/delete_mesage.html"
    model = Inbox
    context_object_name = "msg"
    success_url = reverse_lazy("dashboard")


class UserAddView(generic.CreateView):
    template_name = "pages/add_user.html"
    form_class = UserCreationForm
    success_url = reverse_lazy("dashboard")
    model = User
    context_object_name = "user"


class UserDeleteView(generic.DeleteView):
    template_name = "pages/deleteuser.html"
    success_url = reverse_lazy("dashboard")
    model = User
    context_object_name = "user"


class UserEditView(generic.UpdateView):
    template_name = "accounts/useredit.html"
    form_class = UserEditForm
    success_url = reverse_lazy("dashboard")
    model = User
    context_object_name = "user"


class ProfileView(generic.ListView):
    template_name = "accounts/profile.html"
    model = User
    context_object_name = "user"

    def get_context_data(self, **kwargs):
        context = super(ProfileView, self).get_context_data(**kwargs)
        context['posts'] = Posts.objects.all().filter(author=self.request.user.id)
        return context


class UserUpdateView(generic.UpdateView):
    template_name = "accounts/profile_edit.html"
    form_class = UserProfileEditForm
    model = User

    def get_object(self, queryset=None):
        return self.request.user

    def get_success_url(self):
        pk = self.kwargs["pk"]
        return reverse("profile", kwargs={"pk": pk, "username": self.request.user.username})


class ProfileUpdateView(generic.UpdateView, LoginRequiredMixin):
    template_name = "accounts/profile_edit.html"
    form_class = ProfileEditForm
    model = User

    def get_login_url(self):
        if self.request.user.is_anonymous:
            return reverse_lazy('login')

    def get_object(self, queryset=None):
        return self.request.user.profile

    def get_success_url(self):
        pk = self.kwargs["pk"]
        return reverse("profile", kwargs={"pk": pk, "username": self.request.user.username})


class PasswordChance(PasswordChangeView):
    template_name = "accounts/password.html"
    form_class = PasswordChangeForm

    def get_success_url(self):
        messages.success(self.request, "Şifre değiştirme işlemi başarılı.")
        pk = self.kwargs["pk"]
        return reverse("profile", kwargs={"pk": pk, "username": self.request.user.username})


@login_required
def edit(request, pk, username):
    profil_form = ProfileEditForm(instance=request.user.profile, data=request.POST, files=request.FILES)
    if request.method == 'POST':
        if profil_form.is_valid():
            profil_form.save()
            return reverse("profile", kwargs={"pk": pk, "username": username})
    else:
        profil_form = ProfileEditForm(instance=request.user.profile)

    return render(request, 'accounts/profile_edit.html', {'form': profil_form})
