import os.path
from pathlib import Path

from Bio.Phylo.PAML import codeml

BASE_DIR = Path(__file__).resolve().parent

cml = codeml.Codeml()
cml.alignment = os.path.join(BASE_DIR, "bioinformatic", "files", "aligment.fasta")
cml.tree = os.path.join(BASE_DIR, "bioinformatic", "files", "tree.xml")
cml.ctl_file = os.path.join(BASE_DIR, "bioinformatic", "files", "codeml.ctl")
cml.out_file = os.path.join(BASE_DIR, "bioinformatic", "files", "result_codeml.txt")
cml.working_dir = os.path.join(BASE_DIR, "bioinformatic", "files")
cml_exe = os.path.join(BASE_DIR, "bioinformatic", "apps", "paml4.9j", "bin", "codeml.exe")


results = cml.run(verbose=True, command=cml_exe)
