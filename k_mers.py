from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
# import random
import csv


def make_kmer(filename, formatname, k):
    for seq_record in SeqIO.parse(filename, formatname):
        i = 0
        while True:
            seq = seq_record.seq[i:i+k]
            yield seq
            i += 1


def collect_kmers(filename, formatname, k):
    kmer = make_kmer(filename, formatname, k)
    kmers = set([])
    while True:
        seq = next(kmer)
        if (len(seq) < k or seq.find("N") != -1):
            break
        else:
            kmers.add(seq)
    return list(kmers)


def kmer_index(filename, formatname, kmer):
    for seq_record in SeqIO.parse(filename, formatname):
        return seq_record.seq.find(kmer)


def main():
    # for HOR in [1, 2, 11]:
    #     for k in range(16, 20):
    #         outputf = "chr" + str(HOR) + "_" + str(k) + "merlist.csv"
    #         inputf = "chr" + str(HOR) + "_HOR.fa"
    #         kmer = collect_kmers(inputf, "fasta", k)
    #         with open(outputf, "w") as outfile:
    #             writer = csv.writer(outfile, delimiter=",")
    #             for i in kmer:
    #                 row = ["HOR=" + str(HOR) + ":" + str(kmer_index(inputf, "fasta", i)), i]
    #                 writer.writerow(row)

    cent = [str(x) for x in range(1, 21)]
    cent.extend(["X", "Y"])
    cent.remove("5")
    cent.remove("19")

    for x in cent:
        print(x)
        inputf = "chr" + str(x) + "_HOR.fa"
        kmer = make_kmer(inputf, "fasta", 20)
        while True:
            try:
                test = next(kmer)
            except:
                break
            else:
                print(test)


if __name__ == "__main__":
    main()

# def make_shuffle_seq(record):
#     nuc_list = list(record)
#     random.shuffle(nuc_list)
#     return Seq("".join(nuc_list), record.alphabet)
#
#
# print(kmer)
# l = kmer[0]
# print(l)
# nuc_list = list(l)
# random.shuffle(nuc_list)
# print("".join(nuc_list))
# m = make_shuffle_seq(l)
# print(m)
#
# for seq_record in SeqIO.parse("test2.fa", "fasta"):
#     print(seq_record.seq.count(l))
#     print(seq_record.seq.count(m))
