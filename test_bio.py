from Bio.Seq import Seq

my_seq = Seq("AGTACACTGGT")
print(f"Sequence: {my_seq}")
print(f"Reverse Complement: {my_seq.reverse_complement()}")
