from django.urls import path
from bioinformatic import views

app_name = "bioinformatic"

urlpatterns = [
    path('bioinformatik-labratory-home/', views.BioinformaticHomeView.as_view(), name="home"),
    path('dna-sequence/', views.sekans, name="sekans_analiz"),
    path('translasyon/', views.translation, name="translation"),
    path('download-fasta/', views.fasta_download, name="fasta_download"),
    path('combine-fasta-download/', views.combine_fasta_download, name="combine_fasta_download"),
    path('fasta-writing/', views.fasta_writing, name="fasta_writing"),
    path('fasta-motif/', views.fasta_motif, name="fasta_motif"),
    path('fasta-add/', views.fasta_add, name="fasta_add"),
    path('fasta-combine/', views.fasta_file_combine, name="fasta_file_combine"),
    path('genbank-read/', views.genbank_read, name="genbank_read"),
    path('fasta-file-translate/', views.fasta_file_translate, name="fasta_file_translate"),
    path('fasta-protein-download/', views.fasta_protein_download, name="fasta_protein_download"),
    path('genbank-region/', views.GenBankResultView.as_view(), name="genbank_region"),
    path('genome-diagram/', views.genome_diagram, name="genome_diagram"),
    path('add-enzym/', views.add_extra_enzim, name="add_extra_enzim"),
    path('genbank-region/<int:pk>/<slug:organism>/', views.GenbankDetailView.as_view(), name="genbank_detail"),
    path('delete-genbank/', views.delete_genbank, name="genbank_delete"),
    path('download-genbank/', views.genbank_download, name="genbank_download"),
    path('genbank-writing/', views.genbank_writing, name="genbank_writing"),
    path('local-alignment/', views.local_alignment, name="local_alignment"),
    path('local-alignment-download/', views.local_alignments_download, name="local_alignment_download"),
    path('global-alignment/', views.global_alignment, name="global_alignment"),
    path('multiple-sequence-alignments-analiz/', views.MultipleSeqAlignment, name="multiple_sequence_alignments"),
    path('global-alignment-download/', views.global_alignments_download, name="global_alignment_download"),
    path('fasta-blast-tools/', views.fasta_blast_tools, name="fasta_blast_tools"),
    path('xml-hsp/', views.hsp_download, name="xml_hsp"),
    path('entrez-file-search/', views.entrez, name="entrez_file_search"),
    path('entrez-file-download/', views.entrez_download, name="entrez_file_download"),
    path('pubmed-articles-id/', views.pubmed, name="pubmed_article"),
    path('pubmed-articles-term/', views.medline_article, name="medline_article"),
    path('filogenetik-agac-download/', views.tree_download, name="tree_download"),
    path('swiss-prot/', views.swiss_prot_file, name="swiss_prot_file"),
    path('swiss-prot-url/', views.swiss_prot_url, name="swiss_prot_url"),
    path('swiss_prot-download/', views.swiss_download, name="swiss_prot_download"),
    path('swiss_prot-table/', views.swiss_list_view, name="swiss_prot_list"),
    path('muscle-aligned-download/', views.muscle_aligned_download, name="muscle_aligned_download"),
    path('blast-xml-download/', views.blast_xml_download, name="blast_xml_download"),
    path('blast-hsp-download/', views.blast_hsp_download, name="blast_hsp_download"),
    path('clustal-stats-download/', views.clustal_stats_download, name="clustal_stats_download"),
    path('clustal-scores-download/', views.clustal_scores_download, name="clustal_scores_download"),
    path('motif-download/', views.motif_download, name="motif_download"),
    path('jaspar-motif-download/', views.jaspar_motif_download, name="jaspar_motif_download"),
    path('pssm-download/', views.pssm_download, name="pssm_download"),
    path('pwm-download/', views.pwm_download, name="pwm_download"),
    path('jaspar-motif-create/', views.jaspar_motif_create, name="jaspar_motif_create"),
    path('nucleotid-matrix-positions-download/', views.nucleotid_matrix_positions_download, name="nucleotid_matrix_positions_download"),
    path('clustal-alignment-download/', views.clustal_alignment_download, name="clustal_alignment_download"),
    path('maximum-likelihood/', views.maxlikehood, name="maximum_likelihood"),
    path('maximum-likelihood-download/', views.maximum_likelihood_download, name="maximum_likelihood_download"),
    path('PhyloXML-download/', views.PhyloXML_download, name="PhyloXML_download"),
    path('clustal-omega-alignment/', views.clustalomega_alignment_download, name="clustalomega_alignment_download"),
    path('swiss_prot/<int:pk>/<slug:accessions>/', views.SwissProtDetailView.as_view(), name="swiss_prot_detail"),
]
