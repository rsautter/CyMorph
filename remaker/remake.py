import numpy
import pandas as pd
import sys

l1 = pd.read_csv(sys.argv[1])
l2 = pd.read_csv(sys.argv[2])
for i in range(3,len(sys.argv)):
        ln = pd.read_csv(sys.argv[i])
        l2.append(ln, ignore_index=True)
l2.to_csv("already_done.csv", index=False)
toBeDone = l1[~(l1["dr7objid"].isin(l2["Id"]))]
toBeDone.to_csv("toBeDone.csv",index=False)
#pd.merge(l1,l2,right_on="dr7objid").save_csv("output.csv",index=False)
