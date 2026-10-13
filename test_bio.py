from Bio.Seq import Seq

from Bio import Entrez
Entrez.email = "liyafayazdeen@gmail.com" #i need to tell NCBI who i am to access the database

print("tester: ")
my_seq = Seq("AGTACACTGGT")
print(f"Sequence: {my_seq}")
print(f"Reverse Complement: {my_seq.reverse_complement()}")


stream = Entrez.efetch(db="nucleotide", id="EU490707", rettype="gb", retmode="text")
print(stream.read())


