import csv
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import generic_dna


cent = [str(x) for x in range(1, 21)]

cent.extend(["X", "Y"])

cent.remove("5")
cent.remove("19")


def kmer_notunique(kmer, seq):
    return seq.find(kmer) != -1


def get_uniquekmer(inputseq_rec, *otherseq_recs):
    i = 0
    k = 1
    kmer = inputseq_rec.seq[i:i+k]
    while len(kmer) != 0:
        for otherseq_rec in otherseq_recs:
            while kmer_notunique(kmer, otherseq_rec.seq):
                k = k + 1
                kmer = inputseq_rec.seq[i:i+k]
            else:
                pass
        yield kmer
        i = i + 1
        k = 1
        kmer = inputseq_rec.seq[i:i+k]


# inputf = "chr1_HOR.fa"
# otherf = "chr2_HOR.fa"
# otherf2 = "chr3_HOR.fa"
# for inputseq in SeqIO.parse(inputf, "fasta"):
#     for otherseq in SeqIO.parse(otherf, "fasta"):
#         for otherseq2 in SeqIO.parse(otherf2, "fasta"):
#             kmer = get_uniquekmer(inputseq, otherseq, otherseq2)
#             while True:
#                 print(next(kmer))
#             else:
#                 break

for x in cent:
    print(x)
    othercent = cent[:]
    othercent.remove(x)
    inputf = "chr" + str(x) + "_HOR.fa"
    for inputseq in SeqIO.parse(inputf, "fasta"):
        otherseqs = []
        for y in othercent:
            print(y)
            otherf = "chr" + str(y) + "_HOR.fa"
            for otherseq in SeqIO.parse(otherf, "fasta"):
                otherseqs.append(otherseq)
        print(otherseqs)
        kmer = get_uniquekmer(inputseq, *otherseqs)
        while True:
            try:
                test = next(kmer)
            except:
                break
            else:
                print(test)
