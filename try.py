import os
from pathlib import Path

from Bio.PDB.PDBParser import PDBParser

from Bio.ExPASy import Prosite, ScanProsite

from Bio import ExPASy

from Bio import SwissProt

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'djangoProject\\bioinformatic\\files\\uniprot_sprot.dat')

handle = open(path)

records = SwissProt.parse(handle)


for record in records:
    print(record.gene_name)

