import pandas as pd
import numpy as np


SF1 = [1, 3, 6, 7, 10, 12, 16]
SF2 = [2, 4, 8, 9, 13, 14, 15, 18, 20]
SF3 = [11, 17, "X"]

for k in range(9, 20):
    df = pd.read_csv("centromere" + str(SF1[0]) + "_" + str(k) + "merlist.csv", names=("position", "kmer"))
    for i in range(1, len(SF1)):
        df1 = pd.read_csv("centromere" + str(SF1[i]) + "_" + str(k) + "merlist.csv", names=("position", "kmer"))
        df2 = pd.concat([df, df1])
        df = df2.drop_duplicates(subset="kmer")
    df.to_csv("SF1_centromere_" + str(k) + "merlist.csv")

for k in range(9, 20):
    df = pd.read_csv("centromere" + str(SF2[0]) + "_" + str(k) + "merlist.csv", names=("position", "kmer"))
    for i in range(1, len(SF2)):
        df1 = pd.read_csv("centromere" + str(SF2[i]) + "_" + str(k) + "merlist.csv", names=("position", "kmer"))
        df2 = pd.concat([df, df1])
        df = df2.drop_duplicates(subset="kmer")
    df.to_csv("SF2_centromere_" + str(k) + "merlist.csv")

for k in range(9, 20):
    df = pd.read_csv("centromere" + str(SF3[0]) + "_" + str(k) + "merlist.csv", names=("position", "kmer"))
    for i in range(1, len(SF3)):
        df1 = pd.read_csv("centromere" + str(SF3[i]) + "_" + str(k) + "merlist.csv", names=("position", "kmer"))
        df2 = pd.concat([df, df1])
        df = df2.drop_duplicates(subset="kmer")
    df.to_csv("SF3_centromere_" + str(k) + "merlist.csv")



SFs = [1, 2, 3]

for k in range(9, 20):
    for SF in SFs:
        print(SF)
        otherSF = SFs[:]
        otherSF.remove(SF)
        inputf = "SF" + str(SF) + "_centromere_" + str(k) + "merlist.csv"
        outputf = "SF" + str(SF) + "_cent_original" + str(k) + "merlist.csv"
        df = pd.read_csv(inputf, names=("position", "kmer"))
#    print(df)
        for x in otherSF:
            inputf2 = "SF" + str(x) + "_" + str(k) + "merlist.csv"
            df2 = pd.read_csv(inputf2, names=("position", "kmer"))
            df3 = pd.concat([df, df2])
            mask = df3[df3.duplicated(subset="kmer")]
#        print(len(mask))
            df4 = pd.concat([df, mask])
#        print(df4.duplicated(subset="kmer", keep="last").value_counts()[True])
            df = df4.drop_duplicates(keep=False, subset="kmer")
#    print(df[["kmer"]])
        df[["kmer"]].to_csv(outputf, index=False)
