import csv
from Bio import SeqIO
from Bio import Seq
from k_mers import collect_kmers


chr = [str(x) for x in range(1, 23)]

chr.extend(["X", "Y"])

cent = [str(x) for x in range(1, 21)]

cent.extend(["X", "Y"])

cent.remove("5")
cent.remove("19")


input_q = "test1.fa"
input_s = "test2.fa"

k = 3

kmer = collect_kmers(input_q, "fasta", k)


for i in kmer:
    for seq_record in SeqIO.parse(input_s, "fasta"):
        index = seq_record.seq.find(i)
        split = seq_record.seq.split(i)
        length = [len(j) for j in split]
        position = [sum(length[:(l + 1)]) + l * k for l in range(len(length) - 1)]
        print(i, split, length, index, position)
