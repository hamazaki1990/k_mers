import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from k_mers import collect_kmers
import pandas as pd
import numpy as np

chr = [str(x) for x in range(1, 23)]

chr.extend(["X", "Y"])

cent = [str(x) for x in range(1, 21)]

cent.extend(["X", "Y"])

cent.remove("5")
cent.remove("19")


for i in chr:
    print("chr" + str(i))
    input_s = "chr" + i + "_livingcentmask.fa"
    try:
        genome = next(SeqIO.parse(input_s, "fasta"))
    except:
        break
    else:
        for j in cent:
            print(str(j))
            input_q = "centromere" + str(j) + "_original_kmer.csv"
            output = "chr" + i + "_centromere" + str(j) + "kmer_hit.gff"
            with open(output, "w") as outfile:
                writer = csv.writer(outfile, delimiter=",")
                writer.writerow(["##gff-version 3"])
                df = pd.read_csv(input_q)
                kmer = df.values.tolist()
#                print(kmer)

                for l in kmer:
#                    print(l[0])
                    seq = Seq(l[0], generic_dna)
                    seq_rc = seq.reverse_complement()
                    index = genome.seq.find(seq)
                    index_rc = genome.seq.find(seq_rc)
                    if index == -1 and index_rc == -1:
                        pass
                    elif index_rc == -1:
                        try:
                            HOR = next(SeqIO.parse(input_q, "fasta"))
                        except:
                            break
                        else:
                            split = genome.seq.split(l)
                            length = [len(x) for x in split]
                            position = [sum(length[:(x + 1)]) + x * k for x in range(len(length) - 1)]
                            for m in position:
                                row = ["chr" + str(i), ".", "kmer_hit", m, m + k, ".", "+", str(j)]
                                writer.writerow(row)
                    elif index == -1:
                        try:
                            HOR = next(SeqIO.parse(input_q, "fasta"))
                        except:
                            break
                        else:
                            split = genome.seq.split(seq_rc)
                            length = [len(x) for x in split]
                            position = [sum(length[:(x + 1)]) + x * k for x in range(len(length) - 1)]
                            for m in position:
                                row = ["chr" + str(i), ".", "kmer_hit", m, m + k, ".", "-", str(j)]
                                writer.writerow(row)
                    else:
                        try:
                            HOR = next(SeqIO.parse(input_q, "fasta"))
                        except:
                            break
                        else:
                            split = genome.seq.split(l)
                            length = [len(x) for x in split]
                            position = [sum(length[:(x + 1)]) + x * k for x in range(len(length) - 1)]
                            split_rc = genome.seq.split(seq_rc)
                            length_rc = [len(x) for x in split_rc]
                            position_rc = [sum(length_rc[:(x + 1)]) + x * k for x in range(len(length_rc) - 1)]
                            for m in position:
                                row = ["chr" + str(i), ".", "kmer_hit", m, m + k, ".", "+", str(j)]
                                writer.writerow(row)
                            for m in position_rc:
                                row = ["chr" + str(i), ".", "kmer_hit", m, m + k, ".", "-", str(j)]
                                writer.writerow(row)
