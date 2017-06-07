import pandas as pd
import sys
import numpy

gal1 = pd.read_csv(sys.argv[1])
gal2 = pd.read_csv(sys.argv[2])
pd.merge(gal1, gal2, left_on='Id',right_on='dr7objid').to_csv("merged.csv",index=False)


