from django.contrib import admin
from bioinformatic.models import LabSlideModel, FastaRead, GenbankRead, PubMedArticle, MedlineArticle, \
    SwissProtModel, BigFileUploadModel, MultipleSequenceAlignment, BiologicalResourcesDatabases, \
    FastaDNAMotifModel, RestrictionModel, FileFormat, DiagramModel, GraphicModels
from django.contrib.auth.models import User

# Register your models here.
admin.site.register(MedlineArticle)
admin.site.register(FastaDNAMotifModel)

@admin.register(FileFormat)
class FileFormatModelAdmin(admin.ModelAdmin):
    list_display = ['name']


@admin.register(GraphicModels)
class GraphicModelModelAdmin(admin.ModelAdmin):
    list_display = ['user', 'graph_type', 'format', 'created']
    list_filter = ['user', 'graph_type', 'format', 'created']
    search_fields = ['user', 'graph_type', 'format']


class RestrictionInline(admin.TabularInline):
    model = RestrictionModel


@admin.register(DiagramModel)
class RestrictionAdmin(admin.ModelAdmin):
    inlines = (RestrictionInline,)


@admin.register(BiologicalResourcesDatabases)
class BiologicalResourcesDatabases(admin.ModelAdmin):
    list_display = ['name', 'created']
    search_fields = ['name', 'created']


@admin.register(MultipleSequenceAlignment)
class LabSlideModelAdmin(admin.ModelAdmin):
    list_display = ['user', 'method', 'tree_type', 'molecule_type', 'created']
    list_filter = ['user', 'method', 'tree_type', 'molecule_type', 'created']
    search_fields = ['user', 'method', 'tree_type', 'created']


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
