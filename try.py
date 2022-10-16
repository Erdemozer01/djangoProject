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

records = SeqIO.parse(os.path.join(BASE_DIR, "bioinformatic", "files", "ls_orchid.gbk.txt"), "genbank")

name = []

for record in records:
    name.append(record)

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
    gd_feature_set.add_feature(feature.feature, color=color, label=True)

gd_diagram.draw(
    format="linear",
    orientation="landscape",
    pagesize="A4",
    fragments=4,
    start=0,
    end=len(record),
)

gd_diagram.write("plasmid_linear.pdf", "PDF")

# Import Python modules
from Bio import GenBank
from reportlab.lib import colors
from GenomeDiagram import GDDiagram, GDUtilities

# Load genome annotations from GenBank file
parser = GenBank.FeatureParser()
fhandle = open(file_path, "r")
genbank_entry = parser.parse(fhandle)
fhandle.close()
# Draw linear diagram of CDS features, with GC content graph
gdd = GDDiagram(file_path)
gdt1 = gdd.new_track(4, greytra
CDS’,
scale_fontsize¼3, greytrack_fontsize¼3)
gdt2 ¼ gdd.new_track(6, greytrack¼1, name¼‘Cv
GC
content’,
scale_fontsize¼3, greytrack_fontsize¼3, height¼2)
gdfs ¼ gdt1.new_set(‘feature’)
gdgs ¼ gdt2.new_set(‘graph’)
graphdata1 ¼ GDUtilities.gc_content(genbank_entry.seq, 1000)
graph1 ¼ gdgs.new_graph(graphdata1, ‘GC
content’, style¼‘line’,
colour¼colors.blue, altcolour¼colors.purple)
graph1.linewidth¼1
for feature in genbank_entry.features:
    if feature.type ¼¼ ‘CDS’:
    gdfs.add_feature(feature, colour¼colors.red)
gdd.draw(format¼‘linear’, orientation¼‘landscape’, tracklines¼0,
                                                              pagesize¼‘A6’, fragments¼10, circular¼1)
# Write image as a PNG raster file, and as a PDF vector image
gdd.write(‘example1.png’, ‘PNG’)
gdd.write(‘example1.pdf’, ‘PDF’)
