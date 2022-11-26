import os.path
import subprocess
import sys
from Bio import SeqIO
from Bio.Phylo.PAML import codeml, baseml
from pathlib import Path
from Bio.Phylo.PAML._paml import PamlError

BASE_DIR = Path(__file__).resolve().parent

command = os.path.join(BASE_DIR, "bioinformatic", "apps", "paml4.9j", "bin", "codeml.exe")
ctl_file = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "codeml.ctl")
results_file = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "result.txt")
paml = codeml.Codeml()

paml.alignment = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "admin_out_alignment.fasta")
paml.tree = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "admin_in_file.dnd")
paml.working_dir = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin")
paml.out_file = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "result.txt")
paml.ctl_file = ctl_file

paml.run(command=command, verbose=True)


