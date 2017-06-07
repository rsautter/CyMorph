import numpy
import csv
import sys
import os
import CyMorph
import main as runner
import time
from mpi4py import MPI 
from mpi4py.MPI import ANY_SOURCE
import ConfigParser


def replaceAfterLastDot(fileName, str2Rep):
    splitted = fileName.split('.')
    splitted[len(splitted)-1] = str2Rep
    return ".".join(splitted)

def headNode(comm,rank):
    line = 1
    nprocess = comm.Get_size()
    filesBeingUsed = numpy.array([0 for i in range(nprocess)],dtype=numpy.int32)
    path = "Field/"
    talkingTo = numpy.zeros(1)
    msg = numpy.zeros(1)

    gal = list(csv.reader(open(sys.argv[1], "rb"), delimiter=','))
    ndata = len(gal)
    header = numpy.array(gal[0])
    gal[1:] = sorted(gal[1:],key=lambda l:l[numpy.where(header == "image")[0][0]])


    for line in range(1, ndata):
        #wait someone say something ....
        comm.Recv(talkingTo,source=ANY_SOURCE)
        print("Recieved request from ",talkingTo[0])

        imgIndex = numpy.where(header == "image")[0][0]
        fileName = gal[line][imgIndex].replace(".gz", "")

        # Download the image, if not found
        if not(os.path.isfile(path + fileName)):
            ra = gal[line][numpy.where(header == "ra")[0][0]]
            dec = gal[line][numpy.where(header == "dec")[0][0]]
            run = gal[line][numpy.where(header == "run")[0][0]]
            rerun = gal[line][numpy.where(header == "rerun")[0][0]]
            camcol = gal[line][numpy.where(header == "camcol")[0][0]]
            #field = gal[line][numpy.where(header == "field")[0][0]]
            dr7id = gal[line][numpy.where(header == "dr7objid")[0][0]]
            
            cmd = "wget --inet4-only -r -nd --directory-prefix=Field http://das.sdss.org/raw/"
            cmd += str(run) + "/"
            cmd += str(rerun) + "/corr/"
            cmd += str(camcol) + "/"
            cmd += fileName + ".gz"
            pr = os.popen(cmd)
            print(cmd)
            pr.read()
            # unzip the image
            cmd = "gzip -d " + path + fileName + ".gz"
            pr = os.popen(cmd)
            pr.read()
        print(int(talkingTo[0]),line)
        print( filesBeingUsed)
        filesBeingUsed[int(talkingTo[0])] = line
        #Answering
        msg[0] = line
        comm.Send(msg,dest=talkingTo[0])

    # Job Done! just saying that to other nodes
    for line in range(1, nprocess):
        comm.Recv(talkingTo, source=ANY_SOURCE)
        msg[0] = -1 
        comm.Send(msg,dest=talkingTo[0])
    mergeTables(nprocess)

def mergeTables(nprocess):
    temp = []
    if not(os.path.isfile("output/result.csv")):
        fout=open("output/result.csv","w")
        with open("output/Process0.csv") as f:
            fout.write(f.readline())
    else:
        fout=open("output/result.csv","a")
    for n in range(0,nprocess):
        fileOutput = "output/Process"+str(n)+".csv"
        print("Merging: "+fileOutput)
        f = open(fileOutput,'r')
        f.next()
        for line in f:
            fout.write(line)
        f.close()
        os.remove(fileOutput)
    fout.close()
                
def printHour(t0,t1):
    day = t1.tm_mday-t0.tm_mday
    hour = t1.tm_hour-t0.tm_hour
    minute = t1.tm_min-t0.tm_min
    second = t1.tm_sec-t0.tm_sec
    print("Done in "+str(day)+"d  "+str(hour)+"h "+str(minute)+"m "+str(second)+"s")

def otherNode(comm,rank):
    buff = numpy.zeros(1)
    #reading file list:
    gal = list(csv.reader(open(sys.argv[1], "rb"), delimiter=','))
    header = numpy.array(gal[0])
    gal[1:] = sorted(gal[1:],key=lambda l:l[numpy.where(header == "image")[0][0]])
    path = "Field/"
    saveResult = "output/Process" + str(rank) + ".csv"
    line = 1
  
    while(line > 0):
        #trying to talk with headnode
        buff[0] = rank
        comm.Send(buff,dest=0)
        comm.Recv(buff,source=0)
        line = int(buff[0])
        print("Line", line)
        if(line < 0):
            #Ok, Finishing
            break       
        imgIndex = numpy.where(header == "image")[0][0]
        fileName = gal[line][imgIndex].replace(".gz", "")
        ra = float(gal[line][numpy.where(header == "ra")[0][0]])
        dec = float(gal[line][numpy.where(header == "dec")[0][0]])
        run = gal[line][numpy.where(header == "run")[0][0]]
        rerun = gal[line][numpy.where(header == "rerun")[0][0]]
        camcol = gal[line][numpy.where(header == "camcol")[0][0]]
        #field = gal[line][numpy.where(header == "field")[0][0]]
        dr7id = gal[line][numpy.where(header == "dr7objid")[0][0]]
        try:
            runner.runParallel("ParallelConfig.ini", path, fileName, ra, dec, dr7id,saveResult)
        except Exception as ex:
            with open('errorLog'+str()+'.csv','a') as f:
                f.write(dr7id+','+fileName+'\n')
            print("Got some Error in id "+dr7id+":"+str(ex.args[0]))       

def moleNode(comm,size,rank):
    t0 = time.gmtime()
    gal = list(csv.reader(open(sys.argv[1], "rb"), delimiter=','))
    header = numpy.array(gal[0])
    #gal[1:] = sorted(gal[1:],key=lambda l:l[numpy.where(header == "image")[0][0]])
    path = "Field/"
    saveResult = "output/Process" + str(rank) + ".csv"
    print(header)
    ndata = len(gal)-1
    myiterMin,myiterMax = int(ndata*float(rank-1)/float(size)+ ndata/float(size)),int(ndata*float(rank)/float(size)+ ndata/float(size))
    count = 0
    configFile = ConfigParser.ConfigParser()
    configFile.read('ParallelConfig.ini')
    cleanIt = configFile.getboolean("File_Configuration","cleanit")
    download =  configFile.getboolean("File_Configuration","download")
    
    for i in range(myiterMin,myiterMax):
        print("Rank:",rank," - ", myiterMin,myiterMax)
        print("Line", i)
        line = i+1
        imgIndex = numpy.where(header == "image")[0][0]
        fileName = gal[line][imgIndex].replace(".gz", "")
        ra = float(gal[line][numpy.where(header == "ra")[0][0]])
        dec = float(gal[line][numpy.where(header == "dec")[0][0]])
        run = gal[line][numpy.where(header == "run")[0][0]]
        rerun = gal[line][numpy.where(header == "rerun")[0][0]]
        camcol = gal[line][numpy.where(header == "camcol")[0][0]]
        #field = gal[line][numpy.where(header == "field")[0][0]]
        dr7id = gal[line][numpy.where(map(lambda k:("id" in k) or ("Id" in k), header))[0][0]]
        #gtype = gal[line][numpy.where(header == "Zoo1")[0][0]

        print("Process ",rank,", line:", line, dr7id)

        #download
        if not(os.path.isfile(path + fileName)) and (download):
            cmd = "wget -r -nd --directory-prefix=Field http://das.sdss.org/raw/"
            cmd += str(run) + "/"
            cmd += str(rerun) + "/corr/"
            cmd += str(camcol) + "/"
            cmd += fileName + ".gz"
            print(cmd)
            pr = os.popen(cmd)
            pr.read()
            # unzip the image
            cmd = "gzip -d " + path + fileName + ".gz"
            pr = os.popen(cmd)
            pr.read()
        if not(os.path.isfile(path + fileName)) and not(download):
            continue
        try:
            runner.runParallel("ParallelConfig.ini", path, fileName, ra, dec, dr7id,saveResult)
            print("Done line",line)
        except Exception as ex:
            if not(os.path.exists("errorLog.csv")):
                with open("errorLog.csv", "w") as f:
                    f.write("file,level,msg\n")
                    f.write(str(dr7id)+",indexes,"+str(ex)+"\n")
            else:
                with open("errorLog.csv", "a") as f :
                    f.write(str(dr7id)+",indexes,"+str(ex)+"\n")
            print("Got some Error in id "+dr7id+":"+str(ex.args[0]))
        if(count > 10) and (cleanIt):
            count = 0
            comm.Barrier()
            #cleaning Field folder
            if(rank==0):
                for toRemove in os.listdir("Field/"):
                    os.remove("Field/"+toRemove)
            comm.Barrier()
        count = count + 1
    
    comm.Barrier()
    if(rank==0):
        mergeTables(size)
        t1 = time.gmtime()
        


#@profile
def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # Am I the headnode?
    print("Id:", rank)
    moleNode(comm,size,rank)
   # if rank == 0:
       # yes
   #    headNode(comm,rank)
   # else:
       # no
   #   otherNode(comm,rank)

if __name__ == "__main__":
    main()
