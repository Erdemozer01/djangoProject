from django.contrib import admin
from .models import Profile, UserMessagesModel


@admin.register(Profile)
class ProfileAdmin(admin.ModelAdmin):
    list_display = ['user', 'birth_day', 'created']
    list_filter = ['user']
    search_fields = ['user']
    raw_id_fields = ['user']

    fieldsets = (
        (None, {
            'fields': (
                'user', 'cover', 'avatar', 'phone', 'job', 'about', 'birth_day')
        }),
        ('Sosyal Medya', {
            'classes': ('collapse ', 'extrapretty'),
            'fields': ('facebook', 'twitter', 'instagram',),
        }),
    )


@admin.register(UserMessagesModel)
class UserMessageAdmin(admin.ModelAdmin):
    list_display = ['sender', 'receiver', 'title', 'created']
    list_filter = ['sender']
    search_fields = ['sender']
