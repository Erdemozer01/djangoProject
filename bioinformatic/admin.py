from django.contrib import admin
from bioinformatic.models import Slide, FastaRead, GenbankRead, BlastQuery, BlastHSP

# Register your models here.
admin.site.register(BlastQuery)
admin.site.register(BlastHSP)


@admin.register(Slide)
class Slide(admin.ModelAdmin):
    list_display = ['title']


@admin.register(FastaRead)
class FastaAdmin(admin.ModelAdmin):
    list_display = ['gene']


@admin.register(GenbankRead)
class GenbankAdmin(admin.ModelAdmin):
    list_display = ['gene']