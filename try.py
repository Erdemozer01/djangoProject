import math
import os
from pathlib import Path
from Bio.Phylo import PhyloXML, PhyloXMLIO
from xml.etree import ElementTree
from Bio import Phylo
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent

in_file = os.path.join(BASE_DIR, "media", "MultipleSequenceAlignment", "admin", "admin_tree.xml")

xml_file_read = PhyloXMLIO.parse(in_file)

names = []
branch_length = []
cols = ["name", "branch_length"]
rows = []

for i in xml_file_read:
    for j in i.clade.find_clades():
        rows.append({
            'name': j.name,
            'branch_length': j.branch_length
        })

df = pd.DataFrame(rows, columns=cols)
df.to_csv('output.csv')

for i in range(94):
    print(i)