import numpy
import pandas as pd
import sys

l1 = pd.read_csv(sys.argv[1],header=None)
l2 = pd.read_csv(sys.argv[2])


print(l1)
pd.merge(l1,l2,right_on="dr7objid").save_csv("output.csv",index=False)
