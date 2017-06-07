
import scipy.stats as stats
from GPA import *
from scipy.ndimage.filters import convolve
from scipy.ndimage.filters import median_filter
import numpy
import gridAlg
import ellipse
import indexes
import galaxyIO
import sys
import os
import ConfigParser

cimport numpy
from libc.math cimport sqrt, pow

verbose = 1

def printIfVerbose(string):
    if (verbose > 0):
        print(string)


cdef class CyMorph:
    cdef int Spirality_Ny, Spirality_Nx, Entropy_KFolds, Concentration_Density
    cdef float Ga_Tolerance, Ga_Angular_Tolerance, Ga_Position_Tolerance, sky, smoothingDegradation, butterOrder,minCut,stampSize
    cdef float d1,d2
    cdef int[:] indexes2Evaluate 
    cdef int verbose
    cdef int onlySegmentation
    cdef int errorVar

    def __init__(self):
        self.Spirality_Ny, self.Spirality_Nx = 200, 400
        self.Entropy_KFolds = 150
        self.Ga_Tolerance = 0.03
        self.Ga_Angular_Tolerance = 0.03
        self.Ga_Position_Tolerance =  1.0
        self.Concentration_Density = 100
        self.verbose = True
        self.sky = -1.0
        self.smoothingDegradation = 0.9
        self.butterOrder = 3
        self.minCut = 0.09
        self.d1,self.d2 = -1.0, -1.0
        self.stampSize = 5.0
        self.onlySegmentation = 0
        self.errorVar = 0
        self.indexes2Evaluate = numpy.array([0,0,0,0,0,0,0,0,0,0,0], dtype=numpy.int32)
        
    
    # Sets
    def setSky(self, float value):
        self.sky = value    
    def setGa_Tolerance(self, float value):
        self.Ga_Tolerance = value 
    def setGa_Position_Tolerance(self, float value):
        self.Ga_Position_Tolerance = value 
    def setSpirality_Ny(self, int value):
        self.Spirality_Ny = value 
    def setSpirality_Nx(self, int value):
        self.Spirality_Nx = value
    def setEntropy_KFolds(self, int value):
        self.Entropy_KFolds = value
    def setGa_Angular_Tolerance(self, float value):
        self.Ga_Angular_Tolerance = value
    def setConcentration_Density(self, int value):
        self.Concentration_Density = value
    def setSmoothnessOrder(self, float value):
        self.butterOrder = value
    def setSmoothnessDegree(self,float value):
        self.smoothingDegradation = 1.0-value
    def setConcentrationDists(self,float d1, float d2):
        self.d1 = d1
        self.d2 = d2
    def setOnlySegmentation(self,int value):
        self.onlySegmentation = value
    def setStampSize(self, float value):
        self.stampSize = value

    def setOnlyIndex(self,c):
        self.indexes2Evaluate = numpy.array([0,0,0,0,0,0], dtype=numpy.int32)
        self.setIndexV(c,int(1))

    def isSetIndex(self,indx):
        if(indx=='C') and (self.indexes2Evaluate[0] == 1):
           return True
        elif(indx=='A')and (self.indexes2Evaluate[1] == 1):
           return True
        elif(indx=='S')and (self.indexes2Evaluate[2] == 1):
           return True
        elif(indx=='H')and (self.indexes2Evaluate[3] == 1):
           return True
        elif(indx=='Sp')and (self.indexes2Evaluate[4] == 1):
           return True
        elif(indx=='Ga')and (self.indexes2Evaluate[5] == 1):
           return True
        elif(indx=='C1')and (self.indexes2Evaluate[6] == 1):
           return True
        elif(indx=='C2')and (self.indexes2Evaluate[7] == 1):
           return True
        elif(indx=='CN')and (self.indexes2Evaluate[8] == 1):
           return True
        elif(indx=='A2')and (self.indexes2Evaluate[9] == 1):
           return True
        elif(indx=='A3')and (self.indexes2Evaluate[10] == 1):
           return True
        elif(indx=='S2')and (self.indexes2Evaluate[11] == 1):
           return True
        elif(indx=='S3')and (self.indexes2Evaluate[12] == 1):
           return True
        elif(indx=='OGa') and (self.indexes2Evaluate[13] == 1):
           return True
        else:
           return False

    def setIndexV(self, indx ):
        if(indx=='C'):
           self.indexes2Evaluate[0] = 1
        elif(indx=='A'):
           self.indexes2Evaluate[1] = 1
        elif(indx=='S'):
           self.indexes2Evaluate[2] = 1
        elif(indx=='H'):
           self.indexes2Evaluate[3] = 1
        elif(indx=='Sp'):
           self.indexes2Evaluate[4] = 1
        elif(indx=='Ga'):
           self.indexes2Evaluate[5] = 1
        elif(indx=='C1'):
           self.indexes2Evaluate[6] = 1
        elif(indx=='C2'):
           self.indexes2Evaluate[7] = 1
        elif(indx=='CN'):
           self.indexes2Evaluate[8] = 1
        elif(indx=='A2'):
           self.indexes2Evaluate[9] = 1
        elif(indx=='A3'):
           self.indexes2Evaluate[10] = 1
        elif(indx=='S2'):
           self.indexes2Evaluate[11] = 1
        elif(indx=='S3'):
           self.indexes2Evaluate[12] = 1
        elif(indx=='OGa'):
           self.indexes2Evaluate[13] = 1
        else:
           raise Exception("Unknown index: "+str(indx))
 


    def _replaceAfterLastDot(self,char* fileName,char* str2Rep):
        splitted = fileName.split('.')
        splitted[len(splitted)-1] = str2Rep
        return ".".join(splitted)

    def clearIt(self,fileName,xtraID):
        #os.remove("cutted/"+str(xtraID)+".fit")
        os.remove(str(xtraID)+".cat")
        os.remove(str(xtraID)+"_seg.fits")
        os.remove(str(xtraID)+"_bcg.fits")   

 
    def _runMaskMaker(self,char* path,char* fileName,char* xtraID,float ra, float dec):
        cuttedFile = str(xtraID)+'.fit'

        configFile = ConfigParser.ConfigParser()
        configFile.read('cfg/paths.ini')
        pythonPath = configFile.get("Path","Python")

        #Change directory, execute maskMaker and get back
        localPath = os.getcwd()
        newPath = os.getcwd()+'/maskMaker'
        os.chdir(newPath)
        cmd = pythonPath+" -W\"ignore\" maskmaker_wcut.py ../"+path+fileName+" "+str(xtraID)+" "+str(ra)+" "+str(dec)+" "+str(self.stampSize)+" >> logMaskMaker.txt"
        pr = os.popen(cmd)
        pr.read()
        os.chdir(localPath)
        with open(str(xtraID)+"_log.txt",'r') as eF:
            self.errorVar = int(eF.read())
        os.remove(str(xtraID)+"_log.txt")
        return cuttedFile

    #@profile
    def run(self,char *fpath,char * image,mask_File='',char *saveResult="",float ra=-1.0,float dec=-1.0,calibratedData=False,char* xtraID='', saveFig=True, clear=False,clip=False):
        cdef:
                #image proprieties 
                int width, heigth, it

                #matrices
                float[:,:] notMasked, mask, segmentation, scaleMatrix, segmentationMask, removedGalaxies, zeros
                float[:,:] matInverse, matSmoth, matSmoothness, transformed, transformedEll, transformed2, transformedEll2

                #indexes:
                float a2, a3, s2, s3,h, sp2, sp3, ga, c1, c2
                float sa2, sa3, ss2, ss3,sh, ssp2, ssp3
                float oa2, oa3, os2, os3,oh, osp2, osp3


        results = [xtraID]
        labels = ["Id"]
        printIfVerbose("Running File: "+fpath+image)
        maskFile = mask_File
        fileName = image
        path= fpath
        if (ra != -1.0) and (dec != -1.0) and (clip==True):
            fileName = self._runMaskMaker(path, image,xtraID,ra,dec)
            path = 'cutted/'

        printIfVerbose("Reading File: "+path+fileName)
        notMasked = galaxyIO.readFITSIMG(path+fileName)
        heigth, width = len(notMasked), len(notMasked[0])
        galaxyIO.runSextractor(path,fileName,xtraID)

        if (maskFile == ''):
            printIfVerbose("No mask... Considering every point in the image")
            mask = numpy.array([[0.0 for j in range(width)] for i in range(heigth)],dtype=numpy.float32)
        else:
            printIfVerbose("Reading mask file "+path+maskFile)
            mask = galaxyIO.readFITSIMG(path+maskFile)

        # primeira rodada do sextractor (com a imagem contaminada)
        printIfVerbose("Running Sextractor")
        segmentation = galaxyIO.readFITSIMG(str(xtraID)+"_seg.fits")
        bcg = galaxyIO.readFITSIMG(str(xtraID)+"_bcg.fits")
        dic, data = galaxyIO.readSextractorOutput(str(xtraID)+".cat")
        dicFiltered, dataFiltered = galaxyIO.filterSextractorData(mask, dic, data)
        if clear:
            self.clearIt(fileName,xtraID)
        printIfVerbose("Interpolating ellipse")
        e = ellipse.ellipse(dic,dataFiltered,calibratedData,self.sky)
        segmentationMask, idGalaxy, segMax = gridAlg.filterSegmentedMask(segmentation,e)
        removedGalaxies = gridAlg.removeOtherGalaxies(notMasked, segmentation, idGalaxy)
        newMat, holes = gridAlg.interpolateEllipse(removedGalaxies,e)

        galaxyIO.plotFITS(newMat,"cutted/"+str(xtraID)+".fit")
        printIfVerbose("Running Sextractor again")

        #segunda rodada do sextractor (com a imagem limpa)
        #calibrando pelo background
        converged = False
        bcgW = 32
        galaxyIO.runSextractor("cutted/", str(xtraID)+".fit", xtraID,["BACK_SIZE"],[bcgW])
        segmentation = galaxyIO.readFITSIMG(str(xtraID)+"_seg.fits")
        bcg = galaxyIO.readFITSIMG(str(xtraID)+"_bcg.fits")
        dic, data = galaxyIO.readSextractorOutput(str(xtraID)+".cat")
        dicFiltered, dataFiltered = galaxyIO.filterSextractorData(mask, dic, data)
        if clear:
            self.clearIt(fileName,xtraID)
        if saveFig:
            galaxyIO.plotFITS(bcg,"imgs/bcg"+str(bcgW)+".fit")
        printIfVerbose("Calculating flux profile")
        try:
            dists, concSeq = gridAlg.getAccumulatedLuminosityM(newMat, bcg, e,self.Concentration_Density)
            total, cutDist = indexes.getTotalLum(dists,concSeq,1.0,self.minCut)
        except Exception:
            self.errorVar = 1

        e = ellipse.ellipse(dic,dataFiltered,calibratedData,self.sky)
        segmentationMask, idGalaxy, segMax = gridAlg.filterSegmentedMask(segmentation,e)
        printIfVerbose("Starting enhancing")
        noBCG = numpy.array([[ newMat[i][j]-bcg[i][j] for j in range(len(newMat[i]))] for i in range(len(newMat))], dtype=numpy.float32)
        dx,dy = convolve(noBCG,numpy.array([[-1,-2,-1],\
                                            [ 0, 0, 0],\
                                            [ 1, 2, 1]])),\
                convolve(noBCG,numpy.array([[-1,0,1],\
                                            [-2,0,2],\
                                            [-1,0,1]])) 
        dp,di = convolve(noBCG,numpy.array([[-2,-1,0],\
                                            [-1, 0,1],\
                                            [ 0, 1,2]])),\
                convolve(noBCG,numpy.array([[ 0, 1, 2],\
                                            [-1, 0, 1],\
                                            [-2,-1, 0]])) 
        gradMod = numpy.array([[sqrt(pow(dp[i][j],2.0)+pow(di[i][j],2.0)+pow(dx[i][j],2.0)+pow(dy[i][j],2.0)) for j in range(len(dx[i]))] for i in range(len(dx))],dtype=numpy.float32)
        maxGrad = numpy.max(gradMod)
        gradModF = numpy.array([[ noBCG[i,j]*(gradMod[i,j])/maxGrad for j in range(len(dx[i]))] for i in range(len(dx))], dtype=numpy.float32)
        
        if saveFig:
             galaxyIO.plotFITS(noBCG,"imgs/noBCG.fit")
             galaxyIO.plotFITS(segmentationMask,"imgs/mask.fit")
             galaxyIO.plotFITS(gradModF,"imgs/gradientMod.fit")
        

        # mascaras, escalas, ...
        #scaleMatrix = numpy.array([[e.findScale(float(j),float(i)) for j in range(width)] for i in range(heigth)], dtype=numpy.float32)
        #scaleMatrixN = numpy.array([[scaleMatrix[j, i]/numpy.max(scaleMatrix) for j in range(width)] for i in range(heigth)], dtype=numpy.float32)
        printIfVerbose("Making mask")
        segmentationMask, idGalaxy, segMax = gridAlg.filterSegmentedMask(segmentation,e)
        mat = gridAlg.applyMask(notMasked, mask)
        matSexSeg = gridAlg.applyMask(notMasked, segmentationMask)
        #if(self.Otsu ==1):
        #    otsuMask, otsu = gridAlg.filterOtsu(notMasked, e, self.Entropy_KFolds)
        #    matOtsuSeg = gridAlg.applyMask(notMasked, otsuMask) 

        #making smoothed image
        printIfVerbose("Starting Indexes")
        if(self.isSetIndex("S") or self.isSetIndex("S2") or self.isSetIndex("S3")or self.isSetIndex("S1")):
            printIfVerbose("Smoothing")
            matSmoothness = gridAlg.filterButterworth2D(gradModF,self.smoothingDegradation, self.butterOrder,e)
            diffMat = numpy.array([[(matSmoothness[i][j]-gradModF[i][j]) for j in range(len(newMat[i]))] for i in range(len(newMat))], dtype=numpy.float32)
            maximo = numpy.max(diffMat)
            diffSmoothed = numpy.array([[(matSmoothness[i][j]-noBCG[i][j])/maximo for j in range(len(newMat[i]))] for i in range(len(newMat))], dtype=numpy.float32)
            rotated = indexes.rotateImage(newMat,180)
            diffRotated = numpy.array([[ (newMat[i][j]-rotated[i][j]) for j in range(len(newMat[i]))] for i in range(len(newMat))], dtype=numpy.float32)
            
            if(saveFig):
                galaxyIO.plotFITS(matSmoothness,"imgs/smoothed.fits")
                galaxyIO.plotFITS(diffSmoothed,"imgs/smoothDiff.fits")
                galaxyIO.plotFITS(diffRotated,"imgs/rotateDiff.fits")
            printIfVerbose("Smoothed")   
        printIfVerbose("Object position: "+str((e.posx,e.posy)))

        
        if(self.onlySegmentation == 0):
            printIfVerbose("Without masking")
            if(self.isSetIndex("A") or self.isSetIndex("A2")):
                    a2, matInverse, a2Corr = indexes.asymmetryFunction(gradModF, mask, stats.pearsonr, e)
                    if(saveFig):
                            numpy.savetxt("imgs/asymmetry.txt",numpy.array(a2Corr).T)
                    labels.append("A2")
                    results.append(a2)
            if(self.isSetIndex("A") or self.isSetIndex("A3")):
                    a3, matInverse, a3Corr = indexes.asymmetryFunction(gradModF, mask, stats.spearmanr, e)
                    if(saveFig):
                            numpy.savetxt("imgs/asymmetry.txt",numpy.array(a3Corr).T)
                    labels.append("A3")
                    results.append(a3)
            if(self.isSetIndex("S") or self.isSetIndex("S1")):
                    s1= indexes.smoothness1(gradModF,matSmoothness, mask, e)
                    labels.append("S1")
                    results.append(s1)
            if(self.isSetIndex("S") or self.isSetIndex("S2")):
                    s2, s22RpCorr = indexes.smoothnessFunction(gradModF,matSmoothness, mask, stats.pearsonr,e)
                    if(saveFig):
                            numpy.savetxt("imgs/smoothness.txt",numpy.array(s22RpCorr).T)      
                    labels.append("S2")
                    results.append(s2)
            if(self.isSetIndex("S") or self.isSetIndex("S3")):
                    s3, s32RpCorr = indexes.smoothnessFunction(gradModF,matSmoothness, mask, stats.spearmanr,e)
                    if(saveFig):
                            numpy.savetxt("imgs/smoothness.txt",numpy.array(s32RpCorr).T)  
                    labels.append("S3")
                    results.append(s3)
            if(self.isSetIndex("H")):
                    h = indexes.entropyFunction(mat, self.Entropy_KFolds)
                    if (h>1.0) or (h < 0.0):
                            raise Exception("Unexpected Entropy value:"+str(h))
                    labels.append("H")
                    results.append(h)
            if(self.isSetIndex("Ga")):
                    printIfVerbose("Starting GPA")
                    gpaObject = GPA(noBCG)
                    gpaObject.setPosition(e.posx,e.posy)
                    gpaObject.r = numpy.min(numpy.array([e.maxRad, float(width)/2.0,float(heigth)/2.0]))
                    ga = gpaObject.evaluate(mtol=self.Ga_Tolerance, ftol=self.Ga_Angular_Tolerance, ptol=self.Ga_Position_Tolerance,mask=mask)
                    if (ga>2.0) or (ga < 0.0):
                                raise Exception('Unexpected Ga value:'+str(ga))
                    labels.append("Ga")
                    results.append(ga)
                    if(self.isSetIndex("OGa")):
                                labels.append("OGa")
                                results.append(gpaObject.generate_triangulation_points(self.Ga_Tolerance))

        printIfVerbose("Starting Concentartion")
        if(self.isSetIndex("C") or self.isSetIndex("C1")):
                labels.append("C1")
                if(self.errorVar == 1):
                        results.append(numpy.nan)
                else:
                        c1, total = indexes.concentrationFunction(dists, concSeq, 0.8, 0.2, self.minCut)
                        results.append(c1)
                if(saveFig):
                        normalized = [it/total for it in concSeq]
                        numpy.savetxt("imgs/concentrationN.txt",numpy.array([dists,normalized]).T)
        if(self.isSetIndex("C") or self.isSetIndex("C2")):
                labels.append("C2")
                if(self.errorVar == 1):
                        results.append(numpy.nan)
                else:
                        c2, total = indexes.concentrationFunction(dists, concSeq, 0.9, 0.5, self.minCut)
                        results.append(c2)
                if(saveFig):
                        normalized = [it/total for it in concSeq]
                        numpy.savetxt("imgs/concentrationN.txt",numpy.array([dists,normalized]).T)
        if (self.isSetIndex("C") or self.isSetIndex("CN")) and (self.d1>0.0) and (self.d2>0.0) and (self.d1<1.0)and (self.d2<1.0):
                labels.append("CN")
                if(self.errorVar == 1):
                        results.append(numpy.nan)
                else:
                        cn,total = indexes.concentrationFunction(dists, concSeq, self.d1, self.d2, self.minCut)
                        results.append(cn)
        
        printIfVerbose("Sextractor segmentation")
        # Com segmentacao do sextractor:
        if(self.isSetIndex("A") or self.isSetIndex("A2")):
                sa2, matInverse, a2SexCorr = indexes.asymmetryFunction(gradModF, segmentationMask, stats.pearsonr, e)
                if(saveFig):
                        numpy.savetxt("imgs/sasymmetry.txt",numpy.array(a2SexCorr).T)
                labels.append("sA2")
                results.append(sa2)
        if(self.isSetIndex("A") or self.isSetIndex("A3")):
                sa3, matInverse, a3SexCorr = indexes.asymmetryFunction(gradModF, segmentationMask,  stats.spearmanr, e)
                if(saveFig):
                        numpy.savetxt("imgs/sasymmetry.txt",numpy.array(a3SexCorr).T)
                labels.append("sA3")
                results.append(sa3)
        if(self.isSetIndex("S") or self.isSetIndex("S1")):
                ss1= indexes.smoothness1(gradModF,matSmoothness, segmentationMask, e)
                labels.append("sS1")
                results.append(ss1)
        if(self.isSetIndex("S") or self.isSetIndex("S2")):
                ss2, s2SexCorr = indexes.smoothnessFunction(gradModF,matSmoothness, segmentationMask,  stats.pearsonr,e)
                if(saveFig):
                        numpy.savetxt("imgs/ssmoothness.txt",numpy.array(s2SexCorr).T)
                labels.append("sS2")
                results.append(ss2)
        if(self.isSetIndex("S") or self.isSetIndex("S3")):
                ss3, s3SexCorr = indexes.smoothnessFunction(gradModF,matSmoothness, segmentationMask,  stats.spearmanr,e)
                if(saveFig):
                        numpy.savetxt("imgs/ssmoothness.txt",numpy.array(s3SexCorr).T)
                labels.append("sS3")
                results.append(ss3)
        if(self.isSetIndex("H")):
                sh = indexes.entropyFunction(matSexSeg,self.Entropy_KFolds)
                if (sh>1.0) or (sh < 0.0):
                        raise Exception("Unexpected Sextractor Segmentation Entropy value:"+str(sh))
                labels.append("sH")
                results.append(sh)
        if(self.isSetIndex("Ga")):
                printIfVerbose("Starting GPA")
                gpaObject = GPA(newMat)
                gpaObject.setPosition(e.posx,e.posy)
                gpaObject.r = numpy.min(numpy.array([e.maxRad, float(width)/2.0,float(heigth)/2.0]))
                ga = gpaObject.evaluate(mtol=self.Ga_Tolerance, ftol=self.Ga_Angular_Tolerance, ptol=self.Ga_Position_Tolerance,mask=segmentationMask)
                if (ga>2.0) or (ga < 0.0):
                            raise Exception('Unexpected Ga value:'+str(ga))
                labels.append("sGa")
                results.append(ga)
                if(self.isSetIndex("OGa")):
                            labels.append("OGa")
                            results.append(gpaObject.generate_triangulation_points(self.Ga_Tolerance))   
        
        labels.append("Error")
        results.append(self.errorVar)
        if (len(saveResult) != 0):
            outputVector= [results]
            if not(os.path.isfile(saveResult)):
                numpy.savetxt(saveResult, numpy.array([labels]), delimiter=',', fmt="%s")
            with open(saveResult,'a') as f_handle:
                numpy.savetxt(f_handle, numpy.array(outputVector), delimiter=',', fmt="%s")
            printIfVerbose("File "+fileName+" done.")
        else:
            outputVector= [results]
            numpy.savetxt("Result.csv", numpy.array([labels]), delimiter=',', fmt="%s")
            with open("Result.csv",'a') as f_handle:
                numpy.savetxt(f_handle, numpy.array(outputVector), delimiter=',', fmt="%s")
            printIfVerbose("Results saved in Result.csv")
    
