from Bio import SeqIO
from Bio.SeqUtils import GC

chr = [str(x) for x in range(1, 23)]

chr.extend(["X", "Y"])
GCcont = []

cent = [str(x) for x in range(1, 21)]

cent.extend(["X", "Y"])

cent.remove("5")
cent.remove("19")

for i in cent:
    inputfile = "chr" + str(i) + "_centromere.fa"
    try:
        seq = next(SeqIO.parse(inputfile, "fasta"))
    except:
        break
    else:
        print(inputfile, GC(seq.seq))
        GCcont.append(GC(seq.seq))

for i in chr:
    inputfile = "chr" + str(i) + "_livingcentmask.fa"
    try:
        seq = next(SeqIO.parse(inputfile, "fasta"))
    except:
        break
    else:
        print(inputfile, GC(seq.seq))
        GCcont.append(GC(seq.seq))
print(sum(GCcont)/len(GCcont))

for i in [1, 2, 11]:
    inputfile = "chr" + str(i) + "_HOR.fa"
    try:
        seq = next(SeqIO.parse(inputfile, "fasta"))
    except:
        break
    else:
        print(inputfile, GC(seq.seq))
