import os.path
import subprocess
import sys
from Bio import SeqIO
from Bio.Phylo.PAML import codeml, baseml
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
command = os.path.join(BASE_DIR, "bioinformatic", "apps", "paml4.9j", "bin", "codeml.exe")
ctl_file = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "codeml.ctl")

codeml = codeml.Codeml()

codeml.alignment = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "alignment.fasta")
codeml.tree = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "tree.nwk")
codeml.working_dir = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin")
codeml.out_file = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "result_codeml.txt")
codeml.ctl_file = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "codeml.ctl")

codeml.run(command=command, verbose=True, parse=True)
