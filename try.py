from Bio import Phylo
import os
from pathlib import Path
BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'files\\')

trees = Phylo.read("media/bcl_2.xml", "phyloxml")

Phylo.draw_ascii(trees)
