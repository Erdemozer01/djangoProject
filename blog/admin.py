from django.contrib import admin
from blog.models import Contact, Titles, Cover, Inbox, About, Bottom, Social, Terms, AuthorMessagesModel


# Register your models here.

@admin.register(Titles)
class TitleAdmin(admin.ModelAdmin):
    list_display = ['name', 'type', 'created']


@admin.register(AuthorMessagesModel)
class AuthorMessagesModelAdmin(admin.ModelAdmin):
    list_display = ['sender', 'receiver', 'created']


@admin.register(Terms)
class TermsAdmin(admin.ModelAdmin):
    list_display = ['name', 'created']


@admin.register(Social)
class SocialAdmin(admin.ModelAdmin):
    list_display = ['name', 'created']


@admin.register(About)
class AboutAdmin(admin.ModelAdmin):
    list_display = ['name', 'created']


@admin.register(Contact)
class ContactAdmin(admin.ModelAdmin):
    list_display = ['title', 'email', 'phone', 'address', 'created']


@admin.register(Cover)
class CoverAdmin(admin.ModelAdmin):
    list_display = ['title', 'created']


@admin.register(Bottom)
class BottomAdmin(admin.ModelAdmin):
    list_display = ['title', 'created']


@admin.register(Inbox)
class InboxAdmin(admin.ModelAdmin):
    list_display = ['name', 'email', 'topic', 'phone', 'created']
    list_filter = ['name']
