from Bio import SeqIO
from pathlib import Path
import os
import gzip

BASE_DIR = Path(__file__).resolve().parent.parent
path = os.path.join(BASE_DIR, 'djangoProject\\bioinformatic\\files\\gbvrl2.seq')

reading = open(path, 'r')

records = SeqIO.parse(path, "genbank")

reading.close()

seq = []

for record in records:
    for feature in record.features:
        if feature.qualifiers.get('translation') is not None:
            print(feature.qualifiers.get('protein_id')[0])
            print(feature.qualifiers.get('translation')[0])

print(seq)