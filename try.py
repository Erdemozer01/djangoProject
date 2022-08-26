from Bio.PDB.PDBParser import PDBParser

from Bio.ExPASy import Prosite

from Bio import ExPASy


handle = ExPASy.get_prosite_raw("PS00003")
record = Prosite.read(handle)


print(record.type)
print(record.name)
print(record.created)