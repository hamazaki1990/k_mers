from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.Graphics import BasicChromosome
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio import SeqIO
from BCBio import GFF
import csv


# in_file = "chr9_20mer_chr2HOR_hit_test.gff"
#
# in_handle = open(in_file)
# for rec in GFF.parse(in_handle):
#     print(rec.features[0])
# in_handle.close()

# genome_rec = SeqIO.read("chr2_repeatmasked.fa", "fasta")
# HOR_rec = SeqIO.read("chr2_HOR.fa", "fasta")
#
# with open('chr2_20mer_chr2HOR_hit_test.csv', 'r') as f:
#     reader = csv.reader(f)
#     for row in reader:
#         g_start = ExactPosition(int(row[4]))
#         g_end = ExactPosition(int(row[4]) + 20)
#         H_start = ExactPosition(int(row[1]))
#         H_end = ExactPosition(int(row[1]) + 20)
#         g_feature = SeqFeature(FeatureLocation(g_start, g_end), type="kmer_hit")
#         H_feature = SeqFeature(FeatureLocation(H_start, H_end), type="kmer")
#         genome_rec.features.append(g_feature)
#         HOR_rec.features.append(H_feature)

chr = [str(x) for x in range(1, 23)]

chr.extend(["X", "Y"])

for k in range(15, 21):
# entries = [("Chr 1", "chr1_repeatmasked.fa", "chr1_kmer_test.gff"),
#            ("Chr 2", "chr2_repeatmasked.fa", "chr2_kmer_test.gff"),
#            ("Chr 3", "chr3_repeatmasked.fa", "chr3_kmer_test.gff"),
#            ("Chr 4", "chr4_repeatmasked.fa", "chr4_kmer_test.gff"),
#            ("Chr 5", "chr5_repeatmasked.fa", "chr5_kmer_test.gff"),
#            ("Chr 6", "chr6_repeatmasked.fa", "chr6_kmer_test.gff"),
#            ("Chr 7", "chr7_repeatmasked.fa", "chr7_kmer_test.gff"),
#            ("Chr 8", "chr8_repeatmasked.fa", "chr8_kmer_test.gff"),
#            ("Chr 9", "chr9_repeatmasked.fa", "chr9_kmer_test.gff"),
#            ("Chr 10", "chr10_repeatmasked.fa", "chr10_kmer_test.gff"),
#            ("Chr 11", "chr11_repeatmasked.fa", "chr11_kmer_test.gff"),
#            ("Chr 12", "chr12_repeatmasked.fa", "chr12_kmer_test.gff"),
#            ("Chr 13", "chr13_repeatmasked.fa", "chr13_kmer_test.gff"),
#            ("Chr 14", "chr14_repeatmasked.fa", "chr14_kmer_test.gff"),
#            ("Chr 15", "chr15_repeatmasked.fa", "chr15_kmer_test.gff"),
#            ("Chr 16", "chr16_repeatmasked.fa", "chr16_kmer_test.gff"),
#            ("Chr 17", "chr17_repeatmasked.fa", "chr17_kmer_test.gff"),
#            ("Chr 18", "chr18_repeatmasked.fa", "chr18_kmer_test.gff"),
#            ("Chr 19", "chr19_repeatmasked.fa", "chr19_kmer_test.gff"),
#            ("Chr 20", "chr20_repeatmasked.fa", "chr20_kmer_test.gff"),
#            ("Chr 21", "chr21_repeatmasked.fa", "chr21_kmer_test.gff"),
#            ("Chr 22", "chr22_repeatmasked.fa", "chr22_kmer_test.gff"),
#            ("Chr X", "chrX_repeatmasked.fa", "chrX_kmer_test.gff"),
#            ("Chr Y", "chrY_repeatmasked.fa", "chrY_kmer_test.gff")
#            ]
    for l in [1, 2, 11]:
        entries = [("Chr 1", "chr1_repeatmasked.fa", "chr1_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 2", "chr2_repeatmasked.fa", "chr2_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 3", "chr3_repeatmasked.fa", "chr3_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 4", "chr4_repeatmasked.fa", "chr4_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 5", "chr5_repeatmasked.fa", "chr5_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 6", "chr6_repeatmasked.fa", "chr6_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 7", "chr7_repeatmasked.fa", "chr7_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 8", "chr8_repeatmasked.fa", "chr8_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 9", "chr9_repeatmasked.fa", "chr9_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 10", "chr10_repeatmasked.fa", "chr10_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 11", "chr11_repeatmasked.fa", "chr11_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 12", "chr12_repeatmasked.fa", "chr12_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 13", "chr13_repeatmasked.fa", "chr13_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 14", "chr14_repeatmasked.fa", "chr14_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 15", "chr15_repeatmasked.fa", "chr15_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 16", "chr16_repeatmasked.fa", "chr16_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 17", "chr17_repeatmasked.fa", "chr17_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 18", "chr18_repeatmasked.fa", "chr18_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 19", "chr19_repeatmasked.fa", "chr19_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 20", "chr20_repeatmasked.fa", "chr20_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 21", "chr21_repeatmasked.fa", "chr21_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr 22", "chr22_repeatmasked.fa", "chr22_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr X", "chrX_repeatmasked.fa", "chrX_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"),
                   ("Chr Y", "chrY_repeatmasked.fa", "chrY_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff")
                   ]


        max_len = 248956422
        telomere_length = 1000000
        chr_diagram = BasicChromosome.Organism()
        chr_diagram.page_size = (29.7*cm, 21*cm)

        for index, (name, fastafile, gfffile) in enumerate(entries):
            record = SeqIO.read(fastafile, "fasta")
            length = len(record)
            features = []
            gff_handle = open(gfffile)
            limit_info = dict(
                gff_type=["kmer_hit", "centromere"],)
            for rec in GFF.parse(gff_handle, limit_info=limit_info):
#        print(rec.features[0])
                features = [f for f in rec.features]
#        print([int(f.qualifiers["phase"][0]) for f in features])
            gff_handle.close()
            for f in features:
                if f.type == "centromere":
                    f.qualifiers["color"] = "gainsboro"
                else:
                    f.qualifiers["color"] = str(l) # [int(f.qualifiers["phase"][0])]
            cur_chromosome = BasicChromosome.Chromosome(name)
            cur_chromosome.scale_num = max_len + 2 * telomere_length

            start = BasicChromosome.TelomereSegment()
            start.scale = telomere_length
            cur_chromosome.add(start)

            if features:
                body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
                body.scale = length
                cur_chromosome.add(body)
            else:
                body = BasicChromosome.ChromosomeSegment()
                body.scale = length
                cur_chromosome.add(body)

            end = BasicChromosome.TelomereSegment(inverted=True)
            end.scale = telomere_length
            cur_chromosome.add(end)

            chr_diagram.add(cur_chromosome)

# chr_diagram.draw("HOR_kmer.pdf", "kmer_hit")
            chr_diagram.draw("chr" + str(l) + "HOR_" + str(k) + "mer.pdf", str(k) + "mer_hit")


# gd_diagram = GenomeDiagram.Diagram("chr2")
# gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
# gd_feature_set = gd_track_for_features.new_set()
#
# for feature in genome_rec.features:
#     if feature.type != "kmer_hit":
#         #Exclude this feature
#         continue
#     if len(gd_feature_set) % 2 == 0:
#         color = colors.blue
#     else:
#         color = colors.lightblue
#     gd_feature_set.add_feature(feature, color=color, label=True)
#
# gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4',
#                 fragments=4, start=0, end=len(genome_rec))

# gd_diagram.write("plasmid_linear.png", "PNG")
