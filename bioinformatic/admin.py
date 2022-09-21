from django.contrib import admin
from bioinformatic.models import LabSlideModel, FastaRead, GenbankRead, BlastQuery, PubMedArticle, MedlineArticle, \
    SwissProtModel, BigFileUploadModel, MultipleSequenceAlignment

# Register your models here.
admin.site.register(BlastQuery)
admin.site.register(MedlineArticle)

@admin.register(MultipleSequenceAlignment)
class LabSlideModelAdmin(admin.ModelAdmin):
    list_display = ['user', 'method', 'algoritma', 'created']
    list_filter = ['user', 'method', 'algoritma', 'created']
    search_fields = ['user', 'method', 'algoritma', 'created']

@admin.register(LabSlideModel)
class LabSlideModelAdmin(admin.ModelAdmin):
    list_display = ['title']


@admin.register(FastaRead)
class FastaAdmin(admin.ModelAdmin):
    list_display = ['gene']


@admin.register(GenbankRead)
class GenbankAdmin(admin.ModelAdmin):
    list_display = ['organism', 'protein_id']
    search_fields = ['organism', 'protein_id']


@admin.register(PubMedArticle)
class PubMedArticleAdmin(admin.ModelAdmin):
    list_display = ['link', 'created']
    list_filter = ['email', 'created']


@admin.register(SwissProtModel)
class SwissProtModelAdmin(admin.ModelAdmin):
    list_display = ['accessions', 'organism', 'sequence_length']
    search_fields = ['accessions', 'organism']
    search_help_text = "Erişim Numarası yada organizmaya göre arama yapınız"


@admin.register(BigFileUploadModel)
class BigFileUploadModelAdmin(admin.ModelAdmin):
    list_display = ['big_file', 'created']

