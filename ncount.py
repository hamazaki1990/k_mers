from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
#import random

# my_dna = Seq("ATCGAGAGAAA", generic_dna)
# my_dna
# nuc_list = list(my_dna)
# random.shuffle(nuc_list)
# print(Seq("".join(nuc_list), my_dna.alphabet))

my_dna = Seq("NNAGNNACACNNNNTGGTNN", generic_dna)
my_dna
print(my_dna.count("N"))
print(my_dna.find("N"))
split = my_dna.split("N")
print(split)
length = [len(x) for x in split]
position = [sum(length[:(x + 1)]) + x for x in range(len(length) - 1)]
print(position)

chr = [str(x) for x in range(1, 23)]

chr.extend(["X", "Y"])

for i in chr:
    print("chr" + str(i))
    input_s = "chr" + i + "_peri_5Mbp.fa"
    try:
        peri = next(SeqIO.parse(input_s, "fasta"))
    except:
        break
    else:
        print(peri.seq.count("N")/len(peri.seq))
