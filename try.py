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

print("normalize".title())
print(motif.counts.normalize(pseudocounts=0.4))
print("\n")
pwm = motif.counts.normalize(pseudocounts={"A": 0.6, "C": 0.4, "G": 0.4, "T": 0.6})
print(pwm)
print("pwd:")
print(pwm.consensus)
print("pssm max")
print(motif.pssm.max)
print("pssm min")
print(motif.pssm.min)
print(motif.pssm)
rpssm = motif.pssm
print(rpssm.calculate(Seq("TACACTGCATTACAACCCAAGCATTA")))

excel = pd.DataFrame(motif.pwm).to_excel('pwm.xlsx')

from Bio.motifs import jaspar

jaspar_motif = jaspar.Motif(matrix_id="Erdem", name="Erdem", instances=motif.instances)

jaspar_motif.pseudocounts = jaspar.calculate_pseudocounts(jaspar_motif)

print(jaspar_motif)

print(jaspar_motif.pseudocounts)

from numpy import array

from Bio.Cluster import distancematrix
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO

records = SeqIO.parse(os.path.join(BASE_DIR, "bioinformatic", "files", "opuntia.fasta.txt"), "genbank")

name = []

for record in records:
    name.append(record)

print(name)
record = name[0]

gd_diagram = GenomeDiagram.Diagram(record.description)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

for feature in record.features:
    if feature.type != "gene":
        # Exclude this feature
        continue
    if len(gd_feature_set) % 2 == 0:
        color = colors.blue
    else:
        color = colors.lightblue
    gd_feature_set.add_feature(feature, color=color, label=True)

gd_diagram.draw(
    format="linear",
    orientation="landscape",
    pagesize="A4",
    fragments=4,
    start=0,
    end=len(record),
)
gd_diagram.write("plasmid_circular_nice.pdf", "PDF")
gd_diagram.write("plasmid_linear.png", "png")
from Bio import SeqRecord
for i in name:
    print(i)
    print("\n")