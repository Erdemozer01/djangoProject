import math
import os
from pathlib import Path
from Bio import Phylo, SeqIO
import dash_cytoscape as cyto
from dash import Dash, html, Input, Output
from Bio.Align.Applications import MuscleCommandline
import plotly.graph_objs as go


BASE_DIR = Path(__file__).resolve().parent

in_file = os.path.join(BASE_DIR, "bioinformatic", "files", "ls_orchid.fasta.txt")
out_file = os.path.join(BASE_DIR, "bioinformatic", "files", "aligned.fasta")
scorefile = os.path.join(BASE_DIR, "bioinformatic", "files", "scorefile.txt")
usetree = os.path.join(BASE_DIR, "bioinformatic", "files", "usetree.txt")

muscle_exe = os.path.join(BASE_DIR, "bioinformatic", "apps", "muscle3.8.425_win32.exe")

muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file,
                                 scorefile=scorefile, tree2=usetree)
muscle_cline()

tree = Phylo.read(usetree, "newick")


Phylo.convert(usetree, "newick", "converted_tree.xml", "phyloxml")