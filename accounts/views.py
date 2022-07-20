from blog.models import Inbox
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
from django.conf import settings


class DashBoardView(ListView, LoginRequiredMixin):
    template_name = 'index.html'
    model = Posts
    queryset = Posts.objects.all().order_by("-id")

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(DashBoardView, self).get_context_data(**kwargs)
        context['users'] = User.objects.all()
        context['messages'] = Inbox.objects.all()
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


class ProfileView(generic.DetailView):
    template_name = "accounts/profile.html"
    model = User
    context_object_name = "user"

    def get_context_data(self, **kwargs):
        context = super(ProfileView, self).get_context_data(**kwargs)
        context['object'] = Posts.objects.filter(author_id=self.request.user.id)
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
