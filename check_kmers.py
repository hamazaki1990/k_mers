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

k = 20

for i in chr:
    print("chr" + str(i))
    input_s = "chr" + i + "_repeatmasked.fa"
    output = "chr" + i + "_" + str(k) + "mer_hit.csv"
    try:
        genome = next(SeqIO.parse(input_s, "fasta"))
    except:
        break
    else:
        for j in cent:
            print(str(j))
            input_q = "chr" + j + "_centromere.fa"
            kmer = collect_kmers(input_q, "fasta", k)

            with open(output, "w") as outfile:
                for l in kmer:
                    index = genome.seq.find(l)
                    if index == -1:
                        pass
                    else:
                        row = ["chr" + str(j) + "_centromere", index, "chr" + str(i), str(l), str(genome.seq.count(l)) + "times hit"]
                        writer = csv.writer(outfile, delimiter=",")
                        writer.writerow(row)
