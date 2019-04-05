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
    input_s = "chr" + i + "_livingcentmask.fa"
    output = "chr" + i + "_" + str(k) + "mer_hit.gff"
    try:
        genome = next(SeqIO.parse(input_s, "fasta"))
    except:
        break
    else:
        with open(output, "w") as outfile:
            writer = csv.writer(outfile, delimiter=",")
            writer.writerow(["##gff-version 3"])
            for j in cent:
                print(str(j))
                input_q = "chr" + str(j) + "_centromere.fa"
                kmer = collect_kmers(input_q, "fasta", k)

                for l in kmer:
                    l_rc = l.reverse_complement()
                    index = genome.seq.find(l)
                    index_rc = genome.seq.find(l_rc)
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
                            split = genome.seq.split(l_rc)
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
                            split_rc = genome.seq.split(l_rc)
                            length_rc = [len(x) for x in split_rc]
                            position_rc = [sum(length_rc[:(x + 1)]) + x * k for x in range(len(length_rc) - 1)]
                            for m in position:
                                row = ["chr" + str(i), ".", "kmer_hit", m, m + k, ".", "+", str(j)]
                                writer.writerow(row)
                            for m in position_rc:
                                row = ["chr" + str(i), ".", "kmer_hit", m, m + k, ".", "-", str(j)]
                                writer.writerow(row)
