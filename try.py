import os.path
import subprocess
import sys
from Bio import SeqIO
from Bio.Phylo.PAML import codeml, baseml
from pathlib import Path
from Bio.Phylo.PAML._paml import PamlError
from Bio.PDB import parse_pdb_header
BASE_DIR = Path(__file__).resolve().parent

pdb_file_path = os.path.join(BASE_DIR, "media", "molecule", "admin", "6z2e.cif")
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
parser = PDBParser(PERMISSIVE=1)

mmcif_dict = MMCIF2Dict(pdb_file_path)
print(mmcif_dict)
print(mmcif_dict.keys())
print(mmcif_dict.values())
sc = mmcif_dict["_exptl_crystal.density_percent_sol"]
print(sc)
y_list = mmcif_dict["_atom_site.Cartn_y"]
print(y_list)





