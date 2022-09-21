from django.urls import path
from bioinformatic import views

app_name = "bioinformatic"

urlpatterns = [
    path('home/', views.BioinformaticHomeView.as_view(), name="home"),
    path('dna_seknas/', views.sekans, name="sekans_analiz"),
    path('translasyon/', views.translation, name="translation"),
    path('fasta_read/', views.fasta_read, name="fasta_read"),
    path('lokus_find/', views.GeneRegionView, name="lokus_find"),
    path('delete/', views.delete_fasta, name="delete"),
    path('download_fasta/', views.fasta_download, name="download"),
    path('fasta_writing/', views.fasta_writing, name="fasta_writing"),
    path('fasta_add/', views.fasta_add, name="fasta_add"),
    path('genbank_read/', views.genbank_read, name="genbank_read"),
    path('genbank_region/', views.GenBankResultView.as_view(), name="genbank_region"),
    path('genbank_region/<int:pk>/<slug:organism>/', views.GenbankDetailView.as_view(), name="genbank_detail"),
    path('delete_genbank/', views.delete_genbank, name="genbank_delete"),
    path('download_genbank/', views.genbank_download, name="genbank_download"),
    path('genbank_writing/', views.genbank_writing, name="genbank_writing"),
    path('local_alignment/', views.local_alignment, name="local_alignment"),
    path('local_alignment/download/', views.local_alignments_download, name="local_alignment_download"),
    path('global_alignment/', views.global_alignment, name="global_alignment"),
    path('multiple_sequence_alignments_analiz/', views.MultipleSeqAlignment, name="multiplesequence_alignments"),
    path('multiple_sequence_alignments_sonuçlar/', views.MultipleSequenceAlignmentListView.as_view(),
         name="multiple_sequence_alignments_results"),
    path('global_alignment/download/', views.global_alignments_download, name="global_alignment_download"),
    path('blast/xml_file_read/', views.xml_file, name="xml_file"),
    path('xml_file_result/', views.blast_result_delete, name="xml_file_delete"),
    path('xml_hsp/', views.hsp_download, name="xml_hsp"),
    path('entrez_file_search/', views.entrez, name="entrez_file_search"),
    path('entrez_file_download/', views.entrez_download, name="entrez_file_download"),
    path('pubmed/articles/id/', views.pubmed, name="pubmed_article"),
    path('pubmed/articles/term/', views.medline_article, name="medline_article"),
    path('filogenetik_agac/download/', views.tree_download, name="tree_download"),
    path('swiss_prot/', views.swiss_prot_file, name="swiss_prot_file"),
    path('swiss_prot_url/', views.swiss_prot_url, name="swiss_prot_url"),
    path('swiss_prot_download/', views.swiss_download, name="swiss_prot_download"),
    path('swiss_prot/table/', views.swiss_list_view, name="swiss_prot_list"),
    path('aligned_download/', views.aligned_download, name="aligned_download"),
    path('clustal_stats_download/', views.clustal_stats_download, name="clustal_stats_download"),
    path('clustal_scores_download/', views.clustal_scores_download, name="clustal_scores_download"),
    path('clustal_alignment/', views.clustal_alignment_download, name="clustal_alignment_download"),
    path('swiss_prot/<int:pk><slug:accessions>/', views.SwissProtDetailView.as_view(), name="swiss_prot_detail"),
]
