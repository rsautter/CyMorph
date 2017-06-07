import numpy
import csv
import math
import pandas as pd

gal = list(csv.reader(open("Zoo1_FINAL_PAR.csv","rb"),delimiter=','))
del gal[0]
output = []
#['dr7objid', 'ra', 'dec', 'z', 'petroR50_r', 'deVAB_r', 'seeing_r', 'run', 'camcol', 'rerun', 'field', 'Zoo1']
header = ['dr7objid', 'ra', 'dec', 'z', 'petroR50_r', 'deVAB_r', 'seeing_r', 'run', 'camcol', 'rerun', 'field', 'Zoo1','image']
for obj in gal:
	if (math.pi*math.pow(float(obj[4]),2.0)/float(obj[5]) >=  40.0*math.pi*pow(float(obj[6])/2.0,2.0)):
		output.append(obj)
print("Tamanho da lista para pi * (petroR50_r^2)*(1/deVAB_r)>=40*pi*(f.seeing_r/2)^2: ", len(output))
#print(output)
data = pd.DataFrame(output,columns=header)
data.to_csv("FilteredGalaxyZoo40.csv",index=False)


#few more steps to filter it
spirals = data.where(data["Zoo1"] =='S')
elipticals = data.where(data["Zoo1"] =='E')
spirals.sort("deVAB_r",ascending=False).head(500).to_csv("spirals.csv",index=False)
elipticals.sort("deVAB_r",ascending=False).head(500).to_csv("elipticals.csv",index=False)

