from django.contrib import admin
from .models import Profile

from blog.models.title import Titles

if Titles.objects.exists():
    for name in Titles.objects.all():
        admin.AdminSite.site_header = name.name
else:
    admin.AdminSite.site_header = "Site Yönetimi"


# Register your models here.


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
