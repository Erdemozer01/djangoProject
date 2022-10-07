import os.path
from Bio import SeqIO
from Bio import motifs
from pathlib import Path
from Bio.Seq import Seq

BASE_DIR = Path(__file__).resolve().parent

file_path = os.path.join(BASE_DIR, "bioinformatic", "files", "opuntia.fasta.txt")

reading = SeqIO.parse(file_path, "fasta")

instances = []

for i in reading:
    instances.append(Seq(f"{i.seq[:10]}"))

motif = motifs.create(instances=instances)

motif.name = "Erdem"

print("Motif")
print(motif)

print("Nucleotide matrix position")
print(motif.counts)

print("consensus".title())
print(motif.consensus)

print("anticonsensus".title())
print(motif.anticonsensus)

print("degenerate consensus sequence".title())
print(motif.degenerate_consensus)
import pandas as pd

print("Compute position weight matrices".title())
print(motif.pwm)
print("\n")

print("background".title())
print(motif.mask)
print("\n")

excel = pd.DataFrame(motif.pwm).to_excel('pwm.xlsx')

from Bio.motifs import jaspar

jaspar_motif = jaspar.Motif(matrix_id="Erdem", name="Erdem", instances=motif.instances)

jaspar_motif.pseudocounts = jaspar.calculate_pseudocounts(jaspar_motif)

print(jaspar_motif)

print(jaspar_motif.pseudocounts)