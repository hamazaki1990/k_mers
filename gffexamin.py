from BCBio.GFF import GFFExaminer

for k in range(15, 21):
    for l in [1, 2, 11]:
        entries = ["chr1_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr2_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr3_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr4_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr5_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr6_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr7_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr8_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr9_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr10_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr11_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr12_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr13_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr14_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr15_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr16_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr17_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr18_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr19_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr20_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr21_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chr22_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chrX_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff",
                   "chrY_" + str(k) + "mer_chr" + str(l) + "HOR_hit_test.gff"
                   ]
        for f in entries:
            in_file = f
            examiner = GFFExaminer()
            in_handle = open(in_file)
            print(examiner.available_limits(in_handle))
            in_handle.close()
