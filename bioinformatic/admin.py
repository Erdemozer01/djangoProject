from django.contrib import admin
from bioinformatic.models import Slide, Fasta, Genbank

# Register your models here.
admin.site.register(Fasta)
admin.site.register(Genbank)


@admin.register(Slide)
class Slide(admin.ModelAdmin):
    list_display = ['title']
