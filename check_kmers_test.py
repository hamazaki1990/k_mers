import csv
from Bio import SeqIO
from Bio import Seq


chr = [str(x) for x in range(1, 23)]

chr.extend(["X", "Y"])

k = 10


def make_kmer(filename, formatname, k):
    rec_iter = SeqIO.parse(filename, formatname)
    i = 0
    while True:
        try:
            seq_record = next(rec_iter)
        except StopIteration:
            break
        else:
            seq = seq_record.seq[i:i+k]
            yield seq
            i += 1


for i in chr:
    print(str(i))
    input_q = "human_livingHORs.fa"
    input_s = "chr" + str(i) + "_livingcentmask.fa"

    kmer = make_kmer(input_q, "fasta", k)
    output = "chr" + str(i) + "_" + str(k) + "mer_hit_test.csv"

    for cent in SeqIO.parse(input_s, "fasta"):

        with open(output, "w") as outfile:
            while True:
                try:
                    seq = next(kmer)
                    print(seq)
                except:
                    break
                else:
                    if len(seq) < k:
                        break
                    elif float(seq.count("N")) < 1:
                        index = cent.seq.find(seq)
                        if index == -1:
                            pass
                        else:
                            row = [seq.id, index, str(seq), str(genome.seq.count(l)) + "times hit"]
                            writer = csv.writer(outfile, delimiter=",")
                            writer.writerow(row)
                    else:
                        pass
