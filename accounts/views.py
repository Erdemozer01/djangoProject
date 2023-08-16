from blog.models import Inbox, Terms, Contact, Bottom, About, Cover, Social, Titles
from django.views import generic
from django.urls import reverse_lazy
from django.shortcuts import render, reverse, redirect
from post.models import Posts
from django.contrib.auth.models import User
from django.contrib import messages
from .forms import *
from django.contrib.auth.mixins import LoginRequiredMixin
from .models import Profile, UserMessagesModel
from django.contrib.auth.decorators import login_required
from django.contrib.auth.views import PasswordChangeView
from django.contrib.auth.forms import PasswordChangeForm
from django.views.generic import *
from django.http import HttpResponseRedirect


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
    context_object_name = "messagess"


class MessageDeleteView(generic.DeleteView):
    template_name = "pages/delete_mesage.html"
    model = Inbox
    context_object_name = "msg"
    success_url = reverse_lazy("dashboard")


def dash_user_delete(request, pk):
    if request.user.is_superuser:
        User.objects.get(pk=pk).delete()
        messages.success(request, 'Kullanıcı başarılı bir şekilde silindi')
        return redirect(request.META['HTTP_REFERER'])
    else:
        from django.conf import settings
        messages.error(request, "Lütfen yetkili girişi yapınız.")
        return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))


def settings(request, pk, username):
    user = Profile.objects.get(user=pk)
    profile_form = ProfileEditForm(request.POST or None, request.FILES or None, instance=request.user.profile)
    user_form = UserEditForm(request.POST or None, instance=request.user)
    password_form = PasswordChangeForm(request.user, data=request.POST or None)
    delete_account_form = DeleteAccountForm(request.POST or None)
    if request.method == "POST":
        profile_form = ProfileEditForm(request.POST or None, request.FILES or None, instance=request.user.profile)
        if user_form.is_valid():
            user_form.save()
            messages.success(request, 'Kullanıcı adınız başarılı bir şekilde değişti')
            return redirect(request.META['HTTP_REFERER'])
        elif password_form.is_valid():
            password_form.save()
            from django.contrib.auth import update_session_auth_hash
            update_session_auth_hash(request, password_form.user)
            messages.success(request, 'Şifreniz başarılı bir şekilde değişmiştir')
            return redirect(request.META['HTTP_REFERER'])
        elif profile_form.is_valid():
            profile_form.save()
            messages.success(request, 'Profil bilgileriniz başarılı bir şekilde güncellenmiştir')
            return redirect(request.META['HTTP_REFERER'])
        elif delete_account_form.is_valid():
            confirm = delete_account_form.cleaned_data['confirm']
            if confirm is True:
                user.delete()
                messages.success(request, 'Profil bilgileriniz başarılı silinmiştir')
                return redirect('login')
            else:
                return redirect(request.META['HTTP_REFERER'])
        else:
            profile_form = ProfileEditForm(instance=request.user.profile)
            user_form = UserEditForm(instance=request.user)

    return render(request, 'accounts/settings.html',
                  {'profile_form': profile_form, 'user_form': user_form, 'password_form': password_form,
                   'delete_account_form': delete_account_form})


class UserDeleteView(generic.DeleteView):
    template_name = "pages/deleteuser.html"
    success_url = reverse_lazy("dashboard")
    model = User


class UserEditView(generic.UpdateView):
    template_name = "accounts/useredit.html"
    form_class = UserEditForm
    model = User
    context_object_name = "user"

    def form_valid(self, form):
        instance = form.save(commit=False)
        first_name = instance.first_name
        last_name = instance.last_name
        email = instance.email
        form.save()
        Profile.objects.filter(user=self.request.user).update(first_name=first_name, last_name=last_name, email=email)
        return super(UserEditView, self).form_valid(form)

    def get_success_url(self):
        pk = self.kwargs["pk"]
        return reverse("profile", kwargs={"pk": pk, "username": self.request.user.username})

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(UserEditView, self).get_context_data(**kwargs)
        context['bre'] = "Kullanıcı adını değiştir"
        return context


class ProfileView(generic.DetailView):
    template_name = "accounts/profile.html"
    model = User
    context_object_name = "user"

    def get(self, request, *args, **kwargs):
        if request.user.is_anonymous:
            from django.conf import settings
            messages.error(request, "Lütfen Giriş Yapınız")
            return redirect('%s?next=%s' % (settings.LOGIN_URL, request.path))
        return super(ProfileView, self).get(request, *args, **kwargs)

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

    def form_valid(self, form):
        instance = form.save(commit=False)
        first_name = instance.first_name
        last_name = instance.last_name
        email = instance.email
        form.save()
        User.objects.filter(username=self.request.user).update(first_name=first_name, last_name=last_name, email=email)
        return super(ProfileUpdateView, self).form_valid(form)

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


class UserMessagesListView(ListView):
    template_name = "accounts/user_messages_list.html"
    model = UserMessagesModel

    def get_queryset(self):
        return UserMessagesModel.objects.filter(receiver=self.request.user).order_by('-created')


class UserMessagesDetailView(DetailView):
    template_name = "accounts/user_messages_detail.html"
    model = UserMessagesModel

    def get(self, request, *args, **kwargs):
        UserMessagesModel.objects.update(status='Okundu')
        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['bre'] = "Mesaj detay"
        return context


class UserMessagesDeleteView(DeleteView):
    template_name = "accounts/user_messages_detail.html"
    model = UserMessagesModel

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['bre'] = "Mesaj detay"
        return context


def user_messages_delete(request, pk):
    UserMessagesModel.objects.get(pk=pk).delete()
    messages.success(request, 'Mesaj başarılı şekilde silindi')
    return HttpResponseRedirect(reverse('user_messages',
                                        kwargs={'pk': request.user.pk, 'username': request.user.username}))


def user_mesaasges_delete_all_read(request, username):
    UserMessagesModel.objects.filter(receiver=request.user, status="Okundu").delete()
    messages.success(request, 'Okunmuş tüm mesajlarınız silindi')
    return HttpResponseRedirect(reverse('user_messages',
                                        kwargs={'pk': request.user.pk, 'username': request.user.username}))


def user_sent_message(request, pk, username):
    form = UserMessagesForm(request.POST or None)
    if request.method == "POST":
        if form.is_valid():
            title = form.cleaned_data['title']
            mesaj = form.cleaned_data['message']
            sender = request.user
            receiver = User.objects.get(pk=pk, username=username)
            UserMessagesModel.objects.create(title=title, sender=sender, receiver=receiver, message=mesaj)
            messages.success(request, 'Mesajnız gönderilmiştir...')

            return HttpResponseRedirect(reverse('user_messages',
                                                kwargs={'pk': request.user.pk, 'username': request.user.username}))
        else:
            form = UserMessagesForm()

    return render(request, 'accounts/sent_message.html', {'form': form})


def user_reply_message(request, pk, username, user_pk):
    form = UserMessagesReplyForm(request.POST or None)
    object = UserMessagesModel.objects.get(pk=pk)
    receiver = User.objects.get(pk=user_pk)
    if request.method == "POST":
        if form.is_valid():
            sender = request.user
            message = form.cleaned_data['message']
            UserMessagesModel.objects.create(sender=sender, receiver=receiver, title=object.title, message=message)
            messages.success(request, 'Mesajnız gönderilmiştir...')
            return HttpResponseRedirect(reverse('user_messages',
                                                kwargs={'pk': request.user.pk, 'username': request.user.username}))
        else:
            form = UserMessagesForm()

    return render(request, 'accounts/reply_message.html', {'form': form, 'object': object})
