import os.path

from Bio.Blast import NCBIWWW, NCBIXML
from pathlib import Path
from Bio import SeqIO

BASE_DIR = Path(__file__).resolve().parent

input_fasta_path = os.path.join(BASE_DIR, "bioinformatic", "files", "opuntia.fasta.txt")

records = SeqIO.parse(input_fasta_path, format="fasta")

protein = None

protein_file = open("protein.txt", "w")

for record in records:
    protein_file.write(f"{record.description}")
    protein_file.write("\n")
    protein_file.write("Sekans Uzunlugu: " + f"{len(record.translate().seq)}")
    protein_file.write("\n")
    protein_file.write(f"{record.translate().seq}")
    protein_file.write(2*"\n")
