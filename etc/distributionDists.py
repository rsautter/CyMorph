import hellinger
import sys
import pandas as pd
import numpy as np

def kl(d1, d2):
        dd1,dd2 = d1/sum(d1),d2/sum(d2)
        where = (d1!=0.0) & (d2!=0.0) 
	return sum(dd1[where]*np.log(dd1[where]/dd2[where]))
def hell(d1,d2):
	dd1,dd2 = d1/sum(d1),d2/sum(d2)
	return hellinger.hellinger2(dd1,dd2)
def adist2(d1,d2):
	return abs((np.average(d1)-np.average(d2)))/(np.std(d1)+np.std(d2))

def dists(ell, sp, var):
	freqe, binse = np.histogram(ell[var],bins= 15, normed=True)
	freqs, binss = np.histogram(sp[var],bins= 15, normed=True)
        bins = binse+binss
	freqe, binse = np.histogram(ell[var],bins=bins, normed=True)
	freqs, binss = np.histogram(sp[var],bins=bins, normed=True)	

	output = [kl(freqe,freqs),kl(freqs,freqe),hell(freqe,freqs),adist2(ell[var],sp[var])]
	return output
	
	

data = pd.read_csv(sys.argv[1])
data =data.where(data["Error"] != 2).dropna()

ell = data.where(data["Zoo1"]=="E").dropna()
sp = data.where(data["Zoo1"]=="S").dropna()

metrics = []
header = ["KL(E,S)","KL(S,E)","H(E,S)","N(E,S)"]
rows = []

metrics.append(dists(ell,sp,"sGa"))
rows.append("Ga")

metrics.append(dists(ell,sp,"sA2"))
rows.append("A2")

metrics.append(dists(ell,sp,"sA3"))
rows.append("A3")

metrics.append(dists(ell,sp,"sS2"))
rows.append("S2")

metrics.append(dists(ell,sp,"sS3"))
rows.append("S3")

metrics.append(dists(ell,sp,"sH"))
rows.append("sH")

metrics.append(dists(ell,sp,"C1"))
rows.append("C1")

metrics.append(dists(ell,sp,"C2"))
rows.append("C2")

metrics.append(dists(ell,sp,"CN"))
rows.append("C3")

df =pd.DataFrame(metrics, index=rows,columns=header)
df.to_csv("dmetrics.csv")
print df



