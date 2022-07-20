from django.contrib import admin
from .models import Posts, Category


# Register your models here.
@admin.register(Category)
class PostAdmin(admin.ModelAdmin):
    list_display = ['title']
    list_filter = ['title']


@admin.register(Posts)
class PostAdmin(admin.ModelAdmin):
    list_display = ['title', 'author', 'tags']
    list_filter = ['title', 'author', 'tags']
