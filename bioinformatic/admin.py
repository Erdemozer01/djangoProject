from django.contrib import admin
from bioinformatic.models import Slide, FastaRead, GenbankRead, BlastQuery, PubMedArticle, MedlineArticle

# Register your models here.
admin.site.register(BlastQuery)
admin.site.register(MedlineArticle)


@admin.register(Slide)
class Slide(admin.ModelAdmin):
    list_display = ['title']


@admin.register(FastaRead)
class FastaAdmin(admin.ModelAdmin):
    list_display = ['gene']


@admin.register(GenbankRead)
class GenbankAdmin(admin.ModelAdmin):
    list_display = ['gene']


@admin.register(PubMedArticle)
class PubMedArticleAdmin(admin.ModelAdmin):
    list_display = ['link', 'created']

    list_filter = ['email', 'created']
