import os.path

from Bio.Blast import NCBIWWW, NCBIXML
from pathlib import Path
from Bio import SeqIO

BASE_DIR = Path(__file__).resolve().parent

input_fasta_path = os.path.join(BASE_DIR, "bioinformatic", "files", "opuntia.fasta.txt")

record = next(SeqIO.parse(input_fasta_path, format="fasta"))

result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))

save_file = open("my_blast.xml", "w")

save_file.write(result_handle.read())

save_file.close()

result_handle = open("my_blast.xml")

from Bio import SearchIO

blast_qresults = SearchIO.read("my_blast.xml", "blast-xml")

blast_records = NCBIXML.parse(result_handle)

for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print("\n")
            print("e value:", hsp.expect)
            print("\n")
            print(hsp)

from pydot import call_graphviz

from Bio import Phylo

tree = Phylo.parse("my_blast.xml", "phyloxml")

Phylo.draw(tree)