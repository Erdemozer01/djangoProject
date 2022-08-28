from django.contrib import admin
from bioinformatic.models import Slide, FastaRead, GenbankRead, BlastQuery, PubMedArticle, MedlineArticle, \
    SwissProtModel

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


@admin.register(SwissProtModel)
class SwissProtModelAdmin(admin.ModelAdmin):
    list_display = ['accessions', 'taxonomy_id', 'organism', 'sequence_length']
    search_fields = ['accessions', 'organism']
    search_help_text = "Erişim Numarası yada organizmaya göre arama yapınız"
