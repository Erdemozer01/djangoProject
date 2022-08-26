from Bio.PDB.PDBParser import PDBParser

from Bio.ExPASy import Prosite, ScanProsite

from Bio import ExPASy

sequence = "MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFT" \
           "CRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLSECKTGNGKNYRGTMSKTKN"

handle = ScanProsite.scan(seq=sequence)

result = ScanProsite.read(handle)

for i in result:
    print(i)