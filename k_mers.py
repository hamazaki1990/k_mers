from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def make_kmers(file, format, k):
    for seq_record in SeqIO.parse(file, format):
        i = 0
        while True:
            seq = seq_record.seq[i:i+k]
            if len(seq) < k:
                break
            else:
                yield seq
        i += 1


kmer = make_kmers("test.fa", "fasta", 3)
