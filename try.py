import os.path

from Bio.Blast import NCBIWWW, NCBIXML
from pathlib import Path
from Bio import SeqIO

BASE_DIR = Path(__file__).resolve().parent

input_fasta_path = os.path.join(BASE_DIR, "bioinformatic", "files", "turtles.fasta.txt")

record = next(SeqIO.parse(input_fasta_path, format="fasta"))

result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))

save_file = open("my_blast.xml", "w")

save_file.write(result_handle.read())

save_file.close()

result_handle = open("my_blast.xml")

blast_records = NCBIXML.read(result_handle)


for blast_record in blast_records:
    print(blast_record.hits)


