import csv
import random
from Bio import SeqIO
from Bio.Seq import Seq
from k_mers import collect_kmers


def make_shuffle_seq(record):
    nuc_list = list(record)
    random.shuffle(nuc_list)
    return Seq("".join(nuc_list), record.alphabet)


chr = [str(x) for x in range(1, 23)]

chr.extend(["X", "Y"])

cent = [str(x) for x in range(1, 21)]

cent.extend(["X", "Y"])

cent.remove("5")
cent.remove("19")


for k in range(9, 16):
    for j in [1, 2, 11]:
        print(str(j))
        input_q = "chr" + str(j) + "_HOR.fa"
        output = str(k) + "mer_chr" + str(j) + "HOR_randomcount.csv"
        kmer = collect_kmers(input_q, "fasta", k)
        count = []

        for i in chr:
            print("chr" + str(i))
            input_s = "chr" + i + "_livingcentmask.fa"
            try:
                genome = next(SeqIO.parse(input_s, "fasta"))
            except:
                break
            else:
                for l in kmer:
                    kmercount = genome.seq.count(l)
                    r_1 = make_shuffle_seq(l)
                    r1count = genome.seq.count(r_1)
                    r_2 = make_shuffle_seq(l)
                    r2count = genome.seq.count(r_2)
                    r_3 = make_shuffle_seq(l)
                    r3count = genome.seq.count(r_3)
                    r_4 = make_shuffle_seq(l)
                    r4count = genome.seq.count(r_4)
                    r_5 = make_shuffle_seq(l)
                    r5count = genome.seq.count(r_5)
                    r_6 = make_shuffle_seq(l)
                    r6count = genome.seq.count(r_6)
                    r_7 = make_shuffle_seq(l)
                    r7count = genome.seq.count(r_7)
                    r_8 = make_shuffle_seq(l)
                    r8count = genome.seq.count(r_8)
                    r_9 = make_shuffle_seq(l)
                    r9count = genome.seq.count(r_9)
                    r_10 = make_shuffle_seq(l)
                    r10count = genome.seq.count(r_10)
                    count.append([kmercount, r1count, r2count, r3count, r4count, r5count, r6count, r7count, r8count, r9count, r10count])

        with open(output, "w") as outfile:
            writer = csv.writer(outfile, delimiter=",")
            writer.writerow(["kmer", sum([count[i][0] for i in range(len(count))])])
            for j in range(1, 10):
                row = ["random" + str(j), sum([count[i][j] for i in range(len(count))])]
                writer.writerow(row)
