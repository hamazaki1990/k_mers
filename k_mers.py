from Bio import SeqIO


def collect_kmers(file, k):
    for seq_record in SeqIO.parse(file, "fasta"):
        kmers_list = []
        for i in range(len(seq_record)):
            kmer = seq_record.seq[i:i+k]
            if len(kmer) == k:
                kmers_list.append(kmer)
        return kmers_list

print(collect_kmers("test.fa", 3))
