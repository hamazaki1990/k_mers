from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


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


print(collect_kmers("test2.fa", "fasta", 160))
