import pandas as pd
import numpy as np


SFs = [1, 2, 3]

for k in range(9, 20):
    for SF in SFs:
        print(SF)
        otherSF = SFs[:]
        otherSF.remove(SF)
        inputf = "SF" + str(SF) + "_" + str(k) + "merlist.csv"
        outputf = "SF" + str(SF) + "_original" + str(k) + "merlist.csv"
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
