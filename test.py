from Bio import SearchIO

blast_qresult = SearchIO.read("C:/Users/virtu/Desktop/my_blast.xml", "blast-xml")
print(blast_qresult)


blast_hsp = blast_qresult[0][0]

print(blast_hsp)