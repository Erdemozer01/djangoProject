from Bio import Entrez

Entrez.email = "ozer2466@gmail.com"

handle = Entrez.efetch(db="nucleotide", id="EU490707", rettype="gb", retmode="text")

print(handle.read())


info = Entrez.einfo()

record = Entrez.read(info)

print(record["DbList"])
for i in record["DbList"]:
    print(i)

