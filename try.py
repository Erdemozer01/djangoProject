import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from Bio import motifs

BASE_DIR = Path(__file__).resolve().parent
fasta = os.path.join(BASE_DIR, "bioinformatic", "files", "aligment.fasta")
fasta_reading = SeqIO.parse(fasta, "fasta")

instances2 = []

for seq in fasta_reading:
    instances2.append(seq.seq[:60])

print(instances2)
import pandas as pd

motif = motifs.create(instances2)
print("counts".title())
print(motif.counts)
print("instances".title())
print(motif.instances)
print("consensus".title())
print(motif.consensus)
print("anticonsensus".title())
with open("consensus.txt", "w") as cons_file:
    cons_file.write("Consensus: ")
    cons_file.write(str(motif.consensus) + "\n")

with open("consensus.txt", "a") as cons_file:
    cons_file.write("Anticonsensus: ")
    cons_file.write(str(motif.anticonsensus))

with open("consensus.txt", "a") as cons_file:
    cons_file.write("Anticonsensus: ")
    cons_file.write(str(motif.anticonsensus))

print(motif.anticonsensus)
print("position specific scoring matrices".title())
print(motif.pssm)
print("position weight matrices".title())
print(motif.pwm)
print(motif.mask)

df = pd.DataFrame(motif.pwm).to_excel("pwm.xlsx")

df2 = pd.DataFrame(motif.pssm).to_excel("pssm.xlsx")

df3 = pd.DataFrame(motif.counts).to_excel("counts.xlsx")
