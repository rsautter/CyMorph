import numpy as np
import os
import sys
import time
import pandas as pd

def runTest(n,ndata):
	"""
		runTest:
			 run a sample of ndata galaxies* with 'n' process
			 *The format for the data is: test(ndata)/ellipticals.csv and test(ndata)/spirals.csv
	"""
        os.chdir("../")
	process = os.popen("mpirun -np "+str(n)+" PCyMorph.sh test"+str(ndata)+"/ellipticals.csv")
        process.read()
        process = os.popen("mpirun -np "+str(n)+" PCyMorph.sh test"+str(ndata)+"/spirals.csv")
        process.read()
        os.chdir("etc/")

if __name__ == "__main__":
	"""
		Recieve as input the number of processors	
	"""
	if(len(sys.argv)!=3):
		raise Exception("Wrong number of inputs! An example: python speedup.py 3 1")
	maxProc, inc = int(sys.argv[1]), int(sys.argv[2])
        ndata = 100
        output = pd.DataFrame([[0.0 for i in range(3)] for j in range(1, maxProc+1,inc)])
        output.columns = ["Process","Time","Ndata"]
	for i in range(1, maxProc+1,inc):
		t0 = time.time()
		runTest(i,ndata)
		total = time.time()-t0
		output["Process"][i-1] = i
		output["Time"][i-1] = total
		output["Ndata"][i-1] = 2*ndata
		print("Average Time: ",total/(2*ndata))
		output.to_csv("Speedup.csv",index=False)
