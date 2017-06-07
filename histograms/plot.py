import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import hellinger
import sys
import scipy.stats as stats

plt.rc('text', usetex=True)

def kl(d1, d2):
        dd1,dd2 = d1/sum(d1),d2/sum(d2)
        where = (d1!=0.0) & (d2!=0.0) 
	return sum(dd1[where]*np.log(dd1[where]/dd2[where]))
def hell(d1,d2):
	dd1,dd2 = d1/sum(d1),d2/sum(d2)
	return hellinger.hellinger2(dd1,dd2)
def adist(d1,d2):
	return abs(np.average(d1)/np.std(d1)-np.average(d2)/np.std(d2))
def adist2(d1,d2):
	return abs((np.average(d1)-np.average(d2)))/(np.std(d1)+np.std(d2))

def plot(r,b,bins):
	freqr = plt.hist(r,bins, normed=True, histtype='step',color='r')
	freqb = plt.hist(b,bins, normed=True, histtype='step',color='b')
	print("Kullback-Leiber:",kl(freqr[0],freqb[0]))
	print("Hellinger:", hell(freqr[0],freqb[0]))
	print("Average distance:",adist(r,b))
	print("Average distance2:",adist2(r,b))

r = pd.read_csv(sys.argv[1])
b = pd.read_csv(sys.argv[2])

plt.figure(1)



print("===================================")
bins = np.array(np.arange(0.0,2.0,0.1))
plot(r["sGa"],b["sGa"],bins)
plt.xlabel(r"Algebric $G_a$", fontsize=15)
#plt.ylabel("Frequency")
plt.savefig("AGa.png")

print("")

plt.figure(2)
print("Geometric GPA")
bins = np.array(np.arange(1.85,2.0,0.005))
plot(r["sOGa"],b["sOGa"],bins)
plt.xlabel(r"Geometric $G_a$", fontsize=15)
#plt.ylabel("Frequency")
plt.savefig("GGa.png")
print("")
plt.figure(3)
print("Asymmetry")
bins = np.array(np.arange(0.0,1.0,0.05))
plot(r["sA3"],b["sA3"],bins)
plt.xlabel(r"$A_3$", fontsize=15)
#plt.ylabel("Frequency")
plt.savefig("A3.png")
print("")
plt.figure(4)
print("Clumpiness")
bins = np.array(np.arange(0.0,1.0,0.05))
plot(r["sS3"],b["sS3"],bins)
plt.xlabel(r"$S_3$", fontsize=15)
#plt.ylabel("Frequency")
plt.savefig("S3.png")
print("")
plt.figure(5)
print("Entropy")
bins = np.array(np.arange(0.2,0.85,0.02))
plot(r["sH"],b["sH"],bins)
plt.xlabel("H", fontsize=15)
#plt.ylabel("Frequency")
plt.savefig("H.png")
print("")
plt.figure(6)
print("Concentration(n)")
bins = np.array(np.arange(0.2,0.8,0.02))
plot(r.dropna()["CN"],b.dropna()["CN"],bins)
plt.xlabel(r"$C_{{65\%} \over {25\%}}$", fontsize=15)
#plt.ylabel("Frequency")
plt.savefig("CN.png")
print("===================================")

##########################################
#Panel
##########################################
#fig = plt.figure(7,figsize=(12.8,14.4))
fig = plt.figure(7,figsize=(12.8,14.4))
fig.add_subplot(321)
bins = np.array(np.arange(0.0,2.0,0.1))
plot(r["sGa"],b["sGa"],bins)
plt.xlabel(r"$G_{Alg}$", fontsize=22)
fig.add_subplot(322)
bins = np.array(np.arange(1.85,2.0,0.005))
plot(r["sOGa"],b["sOGa"],bins)
plt.xlabel(r"$G_A$", fontsize=22)
fig.add_subplot(323)
bins = np.array(np.arange(0.0,1.0,0.05))
plot(r["sA3"],b["sA3"],bins)
plt.xlabel(r"$A_3$", fontsize=22)
bins = np.array(np.arange(0.0,1.0,0.05))
fig.add_subplot(324)
plot(r["sS3"],b["sS3"],bins)
plt.xlabel(r"$S_3$", fontsize=22)
fig.add_subplot(325)
bins = np.array(np.arange(0.2,0.85,0.02))
plot(r["sH"],b["sH"],bins)
plt.xlabel("H", fontsize=22)
fig.add_subplot(326)
bins = np.array(np.arange(0.2,0.8,0.02))
plot(r.dropna()["CN"],b.dropna()["CN"],bins)
plt.xlabel(r"$C$", fontsize=17)
plt.tight_layout()
plt.savefig("panel.png")
#plt.show()



