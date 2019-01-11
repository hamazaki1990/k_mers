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


for k in range(15, 36):
    for i in chr:
        print("chr" + str(i))
        input_s = "chr" + i + "_repeatmasked.fa"
        try:
            genome = next(SeqIO.parse(input_s, "fasta"))
        except:
            break
        else:
            for j in cent:
                print(str(j))
                input_q = "chr" + j + "_HOR.fa"
                output = "chr" + i + "_" + str(k) + "mer_chr" + str(j) + "HOR_hit_test.csv"
                kmer = collect_kmers(input_q, "fasta", k)

                with open(output, "w") as outfile:
                    for l in kmer:
                        index = genome.seq.find(l)
                        if index == -1:
                            pass
                        else:
                            split = genome.seq.split(l)
                            length = [len(x) for x in split]
                            position = [sum(length[:(x + 1)]) + x * k for x in range(len(length) - 1)]
                            row = [l]
                            row.extend(position)
                            writer = csv.writer(outfile, delimiter=",")
                            writer.writerow(row)
