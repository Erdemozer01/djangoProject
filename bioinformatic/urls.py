from django.urls import path
from bioinformatic import views

app_name = "bioinformatic"

urlpatterns = [
    path('home/', views.BioinformaticHomeView.as_view(), name="home"),
    path('dna_seknas/', views.sekans, name="sekans_analiz"),
    path('translasyon/', views.translation, name="translation"),
    path('fasta_read/', views.fasta_read, name="fasta_read"),
    path('lokus_find/', views.lokus_find, name="lokus_find"),
    path('delete/', views.delete_fasta, name="delete"),
    path('download_fasta/', views.fasta_download, name="download"),
    path('fasta_writing/', views.fasta_writing, name="fasta_writing"),
    path('fasta_add/', views.fasta_add, name="fasta_add"),
    path('genbank_read/', views.genbank_read, name="genbank_read"),
    path('genbank_region/', views.genbank_region_find, name="genbank_region"),
    path('delete_genbank/', views.delete_genbank, name="genbank_delete"),
]
