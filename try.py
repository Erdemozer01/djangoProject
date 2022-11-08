from Bio import SeqIO
import os
from pathlib import Path

BASEDIR = Path(__file__).parent
in_file = os.path.join(BASEDIR, "bioinformatic", "files", "ls_orchid.fasta.txt")

SeqIO.convert(in_file=in_file, in_format="fasta", out_file="converted_genbank.genbank", out_format="genbank",
              molecule_type="DNA")