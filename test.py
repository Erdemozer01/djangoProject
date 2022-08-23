from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
parser = MMCIFParser()

import easygui

file = easygui.fileopenbox()

mmcif_dict = MMCIF2Dict(f"{file}")

print(file)

structure = parser.get_structure("1fat", f"{file}")

print(structure.parent)

print(mmcif_dict["_exptl_crystal.density_percent_sol"])