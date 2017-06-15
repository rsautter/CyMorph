import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import norm
import pandas as pd

#import rpy2.robjects as robject
#from pyper import *

import ConfigParser
import os
import sys

#Rpy version:
#def runGaussianMetric(fileDist1,fileDist2,metric):
#    kernel = open("measureDistributions.r",'r').read()
#    kernel += "\ngaussianMetric(\""+fileDist1+"\", \""+fileDist2+"\", \""+metric+"\")\n"
#    result = robject.r(kernel)
#    return (float(result[0]),float(result[1]),float(result[2]))

#pyper version
#def runGaussianMetric(fileDist1,fileDist2,metric):
#    kernel = open("measureDistributions.r",'r').read()
#    kernel += "\nv = gaussianMetric(\""+fileDist1+"\", \""+fileDist2+"\", \""+metric+"\")\n"
#    kernel += "write(v,'routput.txt')\n"
#    print('writing kernel')
#    with open("kernel.r", "w") as o:
#        o.write(kernel) 
#    print("Reading R path")
#    configFile = ConfigParser.ConfigParser()
#    configFile.read('cfg/paths.ini')
#    rpath = configFile.get("Path","R")
#    print(rpath+" kernel.r")
#    process = os.popen(rpath+" kernel.r")
#    print(process.read())
#    with open("routput.txt", "r") as i:
#        result = i.read().split(' ')
#    v1,v2,v3 = float(result[0]),float(result[1]),float(result[2])
#    print(v1,v2,v3)
#    return (v1,v2,v3)


#python version (discrete suppervised):
#source:
def hellinger1(p, q):
    return norm(np.sqrt(p) - np.sqrt(q)) /  np.sqrt(2)

def kl(d1, d2):
    dd1,dd2 = d1/sum(d1),d2/sum(d2)
    where = (d1!=0.0) & (d2!=0.0) 
    return sum(dd1[where]*np.log(dd1[where]/dd2[where]))

def adist2(d1,d2):
    return abs((np.average(d1)-np.average(d2)))/(np.std(d1)+np.std(d2))

#python version(discrete suppervised):
def runMetric(r1,b1,idx):
    r = pd.read_csv(r1)[idx]
    b = pd.read_csv(b1)[idx]
    freqr, bins1,tr = plt.hist(r,bins=10)
    freqb, bins2,tr = plt.hist(b,bins=10)
    bins = sorted(np.concatenate((bins1,bins2)))
    freqr, bins1,tr = plt.hist(r,bins=bins)
    freqb, bins2,tr= plt.hist(b,bins=bins)
    return [kl(freqr,freqb),kl(freqb,freqr),hellinger1(freqr,freqb),adist2(r,b)]


def optimizeCN(r1,r2,nsamples,dataFile1,dataFile2, nprocess=2):
    if(nprocess<2):
        raise Exception("You must specify nprocess>1 (at least one headnode, and a worker)")
    hellinger =[]
    kolmogorovlist = []
    deltalist = []
    r1L = []
    r2L = []
    for i in range(0,nsamples+1):
        for j in range(0,nsamples+1):
            lr1 = round(r1[0]+float(i)*(r1[1]-r1[0])/float(nsamples),4) 
            lr2 = round(r2[0]+float(j)*(r2[1]-r2[0])/float(nsamples),4) 
            if(lr1 <= lr2):
                continue
            print("Starting CN (r1, r2)", lr1,lr2)
            print("Nprocess:",nprocess)
            parser = ConfigParser.ConfigParser()
            parser.add_section("File_Configuration")
            parser.add_section("Output_Configuration")
            parser.add_section("Indexes_Configuration")
            parser.set("File_Configuration","Indexes","C")
            parser.set("Output_Configuration","Verbose",False)
            parser.set("Output_Configuration","SaveFigure",False)
            parser.set("Indexes_Configuration","Concentration_Distances",str(lr1)+","+str(lr2))
            parser.set("Indexes_Configuration","Concentration_Density",100)
            with open("ParallelConfig.ini","w") as cfgfile:
                parser.write(cfgfile)
            cmd ="mpirun -np "+str(nprocess)+" PCyMorph.sh "+dataFile1
            process = os.popen(cmd)
            process.read()
            process = os.popen("mv output/result.csv output/r1.csv")
            process.read()
            cmd ="mpirun -np "+str(nprocess)+" PCyMorph.sh "+dataFile2
            process = os.popen(cmd)
            process.read()
            process = os.popen("mv output/result.csv output/r2.csv")
            process.read()
            print("Running metrics")
            try:
                ga, kol, delta = runGaussianMetric("output/r1.csv","output/r2.csv","CN")
                hellinger.append(ga)
                deltalist.append(delta)
                kolmogorovlist.append(kol)
                r1L.append(lr1)
                r2L.append(lr2)
                with open("optimize/CNoptimization.csv",'w') as o:
                    o.write("r1,r2,hellinger,kolmogorov,delta\n")
                    np.savetxt(o, np.array([r1L,r2L,hellinger,kolmogorovlist,deltalist]).T, delimiter=',')
                process = os.popen("mkdir output/CNr"+str(lr1)+"r"+str(lr2))
                process.read()
                process = os.popen("mv output/r1.csv output/r2.csv output/CNr"+str(lr1)+"r"+str(lr2)+"/")
                process.read()
            except:
                print("Error in Cn -> ",lr1,lr2)

def optimizeGa(gaTol,gaATol,nsamples,dataFile1,dataFile2, nprocess=2):
    if(nprocess<2):
        raise Exception("You must specify nprocess>1 (at least one headnode, and a worker)")
    hellinger =[]
    kolmogorovlist = []
    deltalist = []
    phase = []
    mod = []
    for i in range(0,nsamples):
        for j in range(0,nsamples):
            gaMTol = round(gaTol[0]+float(i)*(gaTol[1]-gaTol[0])/float(nsamples),3) 
            gaAngTol = round(gaATol[0]+float(j)*(gaATol[1]-gaATol[0])/float(nsamples),3) 
            print("Starting Ga (Phase, Angular)", gaMTol,gaAngTol)
            print("Nprocess:",nprocess)
            parser = ConfigParser.ConfigParser()
            parser.add_section("File_Configuration")
            parser.add_section("Output_Configuration")
            parser.add_section("Indexes_Configuration")
            parser.set("File_Configuration","Indexes","Ga")
            parser.set("Output_Configuration","Verbose",False)
            parser.set("Output_Configuration","SaveFigure",False)
            parser.set("Indexes_Configuration","Ga_Tolerance",gaMTol)
            parser.set("Indexes_Configuration","Ga_Angular_Tolerance",gaAngTol)
            parser.set("Indexes_Configuration","Ga_Position_Tolerance",0.0)
            with open("ParallelConfig.ini","w") as cfgfile:
                parser.write(cfgfile)
            cmd ="mpirun -np "+str(nprocess)+" PCyMorph.sh "+dataFile1
            process = os.popen(cmd)
            process.read()
            process = os.popen("mv output/result.csv output/r1.csv")
            process.read()
            cmd ="mpirun -np "+str(nprocess)+" PCyMorph.sh "+dataFile2
            process = os.popen(cmd)
            process.read()
            process = os.popen("mv output/result.csv output/r2.csv")
            process.read()
            print("Running metrics")
            try:
                ga, kol, delta = runGaussianMetric("output/r1.csv","output/r2.csv","Ga")
                hellinger.append(ga)
                deltalist.append(delta)
                kolmogorovlist.append(kol)
                mod.append(gaMTol)
                phase.append(gaAngTol)
                with open("optimize/gaoptimization.csv",'w') as o:
                    o.write("phaseT,modT,hellinger,kolmogorov,delta\n")
                    np.savetxt(o, np.array([phase,mod,hellinger,kolmogorovlist,deltalist]).T, delimiter=',')
                process = os.popen("mkdir output/Ga"+str(gaMTol)+"_ang"+str(gaAngTol)+"/")
                process.read()
                process = os.popen("mv output/r1.csv output/r2.csv output/Ga"+str(gaMTol)+"_ang"+str(gaAngTol)+"/")
                process.read()
            except:
                print("Error in Ga -> ",gaMTol,gaAngTol)

def runPCymorph(datafile,nprocess):
    cmd = "mpirun -np "+str(nprocess)+" PCyMorph.sh "+datafile
    process = os.popen(cmd)
    out = process.read()
    

def optimizeEntropy(hm,nsamples,dataFile1,dataFile2, nprocess=2):
    if(nprocess<2):
        raise Exception("You must specify nprocess>1 (at least one headnode, and a worker)")

    metrics = []
    for i in range(0,nsamples):
            hv = int(hm[0]+float(i)*(hm[1]-hm[0])/float(nsamples))
            print("Starting H:",hv)
            print("Nprocess:",nprocess)
            parser = ConfigParser.ConfigParser()
            parser.add_section("File_Configuration")
            parser.add_section("Output_Configuration")
            parser.add_section("Indexes_Configuration")
            parser.set("File_Configuration","Indexes","H")

            parser.set("File_Configuration","cleanit",False)
            parser.set("File_Configuration","download",False)

            parser.set("Output_Configuration","Verbose",False)
            parser.set("Output_Configuration","SaveFigure",False)
            parser.set("Indexes_Configuration","Entropy_Bins",int(hv))
            with open("ParallelConfig.ini","w") as cfgfile:
                parser.write(cfgfile)
            runPCymorph(dataFile1,nprocess)
            process = os.popen("mv output/result.csv output/r1.csv")
            out = process.read()  
            runPCymorph(dataFile2,nprocess)
            process = os.popen("mv output/result.csv output/r2.csv")
            process.read()
            print("Running metric")
            try:
                nm = runMetric("output/r1.csv","output/r2.csv","sH")
                nm.insert(0,hv)
                metrics.append(nm)
                df = pd.DataFrame(metrics)
                df.columns = ["H","kl1","kl2","hell","N"]
                df.to_csv("optimize/entropy.csv", index=False)
                print(df)
                process = os.popen("mkdir output/Hbin"+str(hv))
                process.read()
                process = os.popen("mv output/r1.csv output/r2.csv output/Hbin"+str(hv)+"/")
                process.read()
            except:
                print("Error in h -> ",hv)

def optimizeSmoothness(sm,nsamples,dataFile1,dataFile2, nprocess=2):
    metricsS2,metricsS3 = [],[]

    for i in range(0,nsamples):
            cv = round(sm[0]+float(i)*(sm[1]-sm[0])/float(nsamples),3)
            print("Starting C:",cv)
            print("Nprocess:",nprocess)
            parser = ConfigParser.ConfigParser()
            parser.add_section("File_Configuration")
            parser.add_section("Output_Configuration")
            parser.add_section("Indexes_Configuration")
            parser.set("File_Configuration","Indexes","S")
            parser.set("Output_Configuration","Verbose",False)
            parser.set("Output_Configuration","SaveFigure",False)
            parser.set("Indexes_Configuration","smooth_degree",cv)
            parser.set("Indexes_Configuration","butterworth_order",2.0)
            parser.set("File_Configuration","cleanit",False)
            parser.set("File_Configuration","download",False)

            with open("ParallelConfig.ini","w") as cfgfile:
                parser.write(cfgfile)
            runPCymorph(dataFile1,nprocess)
            process = os.popen("mv output/result.csv output/r1.csv")
            out = process.read()  
            runPCymorph(dataFile2,nprocess)
            process = os.popen("mv output/result.csv output/r2.csv")
            process.read()
            print("Running metric")
            try:
                nm = runMetric("output/r1.csv","output/r2.csv","sS2")
                nm.insert(0,cv)
                metricsS2.append(nm)
                df = pd.DataFrame(metricsS2)
                df.columns = ["S2","kl1","kl2","hell","N"]
                print(df)
                df.to_csv("optimize/s2.csv", index=False)
            except:
                print("Error in s2 -> ",cv)
            try:
                nm = runMetric("output/r1.csv","output/r2.csv","sS3")
                nm.insert(0,cv)
                metricsS3.append(nm)
                df = pd.DataFrame(metricsS3)
                df.columns = ["S3","kl1","kl2","hell","N"]
                print(df)
                df.to_csv("optimize/s3.csv", index=False)
                process = os.popen("mkdir output/S"+str(cv))
                process.read()
                process = os.popen("mv output/r1.csv output/r2.csv output/S"+str(cv)+"/")
                process.read()
            except:
                print("Error in s3 ->",cv)
            

##The files must already be in Field/
if __name__ == "__main__":
    n=int(sys.argv[1])
    #sm = [0.1,1.0],nsamples=18
    #optimizeSmoothness(sm = [0.1,0.5],nsamples=5,dataFile1="test100/spirals100.csv",dataFile2="test100/ellipticals100.csv",nprocess=n)
    #optimizeCN(r1 = [0.45,0.95],r2 = [0.05,0.55],nsamples=5,dataFile1="test100/spirals100.csv",dataFile2="test100/ellipticals100.csv",nprocess=n)
    #optimizeGa(gaTol=[0.00,0.02],gaATol=[0.00,0.04],nsamples=4,dataFile1="test100/spirals100.csv",dataFile2="test100/ellipticals100.csv",nprocess=n)
    optimizeEntropy(hm = [100,250],nsamples=15,dataFile1="test100/spirals.csv",dataFile2="test100/ellipticals.csv",nprocess=n)
    #print(runGaussianMetric("output/ellipticals.csv","output/spirals.csv","S2"))
