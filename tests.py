import galaxyIO as gio
import indexes as par
import gridAlg
import numpy
import math
import ellipse
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
import galaxyIO
import CyMorph as morf
from scipy.optimize import curve_fit
from scipy import fftpack
from math import sqrt, exp,sin,cos
from scipy.ndimage.filters import convolve


def scatterPlot(x,y,fig,minv,maxv,color='b'):
    plt.figure(fig,figsize=(8, 8))
    plt.plot(x,y,'ko', markersize=.5,color=color)
    plt.xticks(numpy.arange(minv-.5, maxv+.5,1), rotation='vertical')
    plt.yticks(numpy.arange(minv-.5, maxv+.5,1))
    plt.xlim(minv, maxv)
    plt.ylim(minv, maxv)
    plt.axis('off')
    
    plt.grid(True)

def entropyTestMaxMin(fileSave, kfold=4,mlen=100):
    entropies = []
    probability = []
    cteMat = numpy.array([[0.0 for i in range(mlen)] for j in range(mlen)])
    for i in range(mlen /kfold):
        for j in range(kfold):
            cteMat[:,i+int(float(j)*float(mlen)/float(kfold))]=float(j+1)

    print(len(numpy.where(cteMat==1.0)[0])/float(mlen*mlen))
    for i in range(len(cteMat)/kfold):
        entropies.append(par.entropyFunction(cteMat.astype(numpy.float32),kfold))
        probability.append(len(numpy.where(cteMat==1.0)[0])/float(mlen*mlen))
        for j in range(1,kfold):
            cteMat[:,i+int(j*len(cteMat)/kfold)]=1.0
    entropies.append(par.entropyFunction(cteMat.astype(numpy.float32),kfold))
    probability.append(len(numpy.where(cteMat==1.0)[0])/float(mlen*mlen))
    numpy.savetxt(fileSave,numpy.array([entropies,probability]))
    plt.plot(probability,entropies)
    plt.xlabel("p(x=1)")
    plt.ylabel("H(m)")
    plt.ylim(-0.01,1.01)


def skyBoxTest(img):
    print("Testing box ", len(img)/4.,len(img)/4.)
    gridAlg.measureSky(img,len(img)/4.,len(img)/4.,50)
    print(" ")
    
    print("Testing box ",len(img)/4.,3.*len(img)/4.)
    gridAlg.measureSky(img,len(img)/4.,3.*len(img)/4.,50)
    print(" ")
    
    print("Testing box ",3.*len(img)/4.,len(img)/4.)
    gridAlg.measureSky(img,3.*len(img)/4.,len(img)/4.,50)
    print(" ")
    
    print("Testing box ",3.*len(img)/4.,3.*len(img)/4.)
    gridAlg.measureSky(img,3.*len(img)/4.,3.*len(img)/4.,50)

def samplingErrorTest(fileName,density=50,md=100,nsamples=100):
    ones = numpy.array([[1.0 for x in range(md)] for y in range(md)]).astype(numpy.float32)
    zeros = numpy.array([[0.0 for x in range(md)] for y in range(md)]).astype(numpy.float32)
    e.posx, e.posy = md/2,md/2
    error50 = [0.0 for i in range(len(ones)/2)]

    for i in range(nsamples):
        dists, acc = gridAlg.getAccumulatedLuminosityM(ones,zeros,e,nsamples=density)
        error = [error50[i]+ abs(acc[i]- math.pi*(float(i)**2.0))/max(1.0,math.pi*(float(i)**2.0)) for i in range(len(acc))]
        error50 = error
    error50 = numpy.array(error50)/float(nsamples)
    numpy.savetxt(fileName,error50)

def entropyBinTest(img, folds2Test):
    seq = []
    for f in folds2Test:
        print("Testing fold: "+str(f))
        seq.append(par.entropyFunction(img,f))
    return numpy.array(seq)

def smoothingTest(img):
    
    smoothed90 = gridAlg.filterButterworth2D(img,0.9,1.0,len(img)/2,len(img)/2)
    print("90%")
    smoothed70 = gridAlg.filterButterworth2D(img,0.7,1.0,len(img)/2,len(img)/2)
    print("70%")
    smoothed50 = gridAlg.filterButterworth2D(img,0.5,1.0,len(img)/2,len(img)/2)
    print("50%")
    smoothed30 = gridAlg.filterButterworth2D(img,0.3,1.0,len(img)/2,len(img)/2)
    print("30%")
    smoothed10 = gridAlg.filterButterworth2D(img,0.1,1.0,len(img)/2,len(img)/2)
    print("10%")
    
    gio.plotFITS(smoothed90, "tests/smoothness/c0.9.fits")
    gio.plotFITS(smoothed70, "tests/smoothness/c0.7.fits")
    gio.plotFITS(smoothed50, "tests/smoothness/c0.5.fits")
    gio.plotFITS(smoothed30, "tests/smoothness/c0.3.fits")
    gio.plotFITS(smoothed10, "tests/smoothness/c0.1.fits")
    
    p0, = plt.plot(numpy.array(img)[len(img)/2,:],'b+-', label='Original')
    p1, = plt.plot(numpy.array(smoothed90)[len(img)/2,:],'k', label='c = 0.9')
    p2, = plt.plot(numpy.array(smoothed70)[len(img)/2,:],'r', label='c = 0.7')
    p3, = plt.plot(numpy.array(smoothed50)[len(img)/2,:],'g', label='c = 0.5')
    p4, = plt.plot(numpy.array(smoothed30)[len(img)/2,:],'b', label='c = 0.3')
    p5, = plt.plot(numpy.array(smoothed10)[len(img)/2,:],'m', label='c = 0.1')
    plt.xlabel('Position (pixel)')
    plt.ylabel('Instensity')
    
    plt.legend([p0,p1,p2,p3,p4,p5],['original','c = 0.9','c = 0.7','c = 0.5','c = 0.3','c = 0.1'])
    plt.show()


def concentrationTotalTest(testingImg,e,nsample,foldSize):
    zeros = numpy.array([[1000.0 for x in range(len(testingImg[y]))] for y in range(len(testingImg))]).astype(numpy.float32)
    dists, conc = gridAlg.getAccumulatedLuminosityM(testingImg, zeros, e, nsample)


    s5_Sigma_1_Delta = par.getTotalLum(dists,conc,0.5,0.01)
    s1_Sigma_1_Delta = par.getTotalLum(dists,conc,1.0,0.01)
    s15_Sigma_1_Delta = par.getTotalLum(dists,conc,1.5,0.01)

    s5_Sigma_3_Delta = par.getTotalLum(dists,conc,0.5,0.3)
    s1_Sigma_3_Delta = par.getTotalLum(dists,conc,1.0,0.3)
    s15_Sigma_3_Delta = par.getTotalLum(dists,conc,1.5,0.3)
    
    dconc = numpy.array([0.0 for i in range(len(conc)-2)])
    dists2 =  numpy.array([0.0 for i in range(len(conc)-2)])
    for i in range(1,len(conc)-1):
        dconc[i-1] = (conc[i+1] - conc [i-1]) / (dists[i+1]-dists [i-1])
        print(conc[i+1], conc [i-1], conc[i+1] - conc [i-1], (dists[i+1]-dists [i-1])) 
        dists2[i-1] = dists[i]

    dists = numpy.asarray(dists)
    dists2 = numpy.asarray(dists2)
    conc = numpy.asarray(conc)
    dconc = numpy.asarray(dconc)
    
    foldDists = []
    mean, median, maxi, mini,qmax, qmin = [],[],[],[],[],[]
    for i in range(0,len(dconc),foldSize):
        soma = 0.0
        count = 0
        fold = []
        for j in range(i,min(i+foldSize,len(dconc))):
            fold.append(dconc[j])
            soma += dists[j]
            count =count+1
        print(soma, count)
        soma = soma/float(count) if count>0 else soma 
        fold = numpy.array(fold)
        qqmax = numpy.percentile(fold, 75)
        qqmin = numpy.percentile(fold,25)
        cutted = fold[numpy.where((fold>qqmin-1.5*(qqmax-qqmin)) & (fold<qqmax+1.5*(qqmax-qqmin)))]
        fold =numpy.array(fold)
        mean.append(numpy.average(fold))
        median.append(numpy.median(fold))
        maxi.append(numpy.max(cutted))
        mini.append(numpy.min(cutted))
        qmax.append(qqmax)
        qmin.append(qqmin)
        foldDists.append(soma)
            
   
    numpy.savetxt("tests/concentration/conc.csv",numpy.array([dists,conc]).T,fmt='%.5f')
    numpy.savetxt("tests/concentration/Dconc.csv",numpy.array([dists2,dconc]).T,fmt='%.5f')
    
    box = numpy.array([mean,median,maxi,mini,qmax,qmin,foldDists]).T
    header = "mean median max min qmax qmin dist\n"
    with open("tests/concentration/BoxPlot.csv", 'wb') as f:
        f.write(header)
        numpy.savetxt(f,box,fmt='%.5f')

    #numpy.savetxt("tests/concentration/N.0.5.Sigma.0.01.Delta.csv",numpy.array([dists,conc/s5_Sigma_1_Delta]))
    #numpy.savetxt("tests/concentration/N.1.0.Sigma.0.01.Delta.csv",numpy.array([dists,conc/s1_Sigma_1_Delta]).T,fmt='%.5f')
    #numpy.savetxt("tests/concentration/N.1.5.Sigma.0.01.Delta.csv",numpy.array([dists,conc/s15_Sigma_1_Delta]).T,fmt='%.5f')
    #numpy.savetxt("tests/concentration/N.0.5.Sigma.0.3.Delta.csv",numpy.array([dists,conc/s5_Sigma_3_Delta]).T,fmt='%.5f')
    #numpy.savetxt("tests/concentration/N.1.0.Sigma.0.3.Delta.csv",numpy.array([dists,conc/s1_Sigma_3_Delta]).T,fmt='%.5f')
    #numpy.savetxt("tests/concentration/N.1.5.Sigma.0.3.Delta.csv",numpy.array([dists,conc/s15_Sigma_3_Delta]).T,fmt='%.5f')
    
    
def rotationTest(img):
    print("Rotating image")
    imgRotated = par.rotateImage(img, 180.0)
    imgReflected = numpy.array([i[:] for i in img])
    for y in range(len(img)):
        for x in range(len(img[y])):
            imgReflected[y,x] = img[(len(img)-1-y),(len(img[y])-1-x)]
    gio.plotFITS(imgRotated,"tests/Rotation180.fits")
    gio.plotFITS(imgReflected,"tests/Reflection.fits")
    gio.plotFITS(imgRotated-imgReflected,"tests/Rotated-Reflection.fits")
    print("Rotated")

def saveButterworthProfile(n):
    profile = []
    d0 = 50
    seq = []
    for d in range(0, 100):
        seq.append(d)
        profile.append(gridAlg.butterworth(d,d0,n))
    numpy.savetxt("tests/butterworthProfile"+str(n)+".csv", numpy.array([seq, profile]).T)

def testDistances(mat, ellipse, ):
    color=cm.rainbow(numpy.linspace(0,1,len(distances)))
    tseq = []
    nsamplePerDist=1000
    distances=numpy.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8],dtype=numpy.float32)
    for d in distances:
        seq, trash = gridAlg.getSampleDist(mat,nsamplePerDist,d,ellipse)
        avg = numpy.average(seq)
        fft = fftpack.fft(seq/avg)
        tseq=numpy.concatenate([tseq,numpy.real(fft[2:(len(fft)-1)])])
	plt.plot(fft[2:(len(fft)-1)])	
    return max(numpy.abs(tseq))
    
            
   
if __name__ == "__main__":	
    testingImg = "ell.fit" #"img1000_n6_re10_band_rNoise_sim_0.fits"#
    testingImgSp = "sp.fit"
    #testingImgSp ="5.87722952767439E+017.fit"
    path = "tests/input/spiral/"
    pathSp = "tests/input/spiral/"
    #pathSp = "cutted/"

    # Simulated:
    ellipticaWN = gio.readFITSIMG("tests/input/elliptical_Simulated/img1000_n6_re10_band_rNoise_sim_0.fits").astype(numpy.float32)
    spiralWN = gio.readFITSIMG("tests/input/spiral_Simulated/ima1.fits").astype(numpy.float32)

    
    img = gio.readFITSIMG(path+testingImg).astype(numpy.float32)
    imgSp = gio.readFITSIMG(pathSp+testingImgSp).astype(numpy.float32)

    galaxyIO.runSextractor(path,testingImg,0)
    dic, data = galaxyIO.readSextractorOutput("0.cat")
    dic,data = galaxyIO.filterSextractorData(img, dic,data)
    bcg = gio.readFITSIMG("0_bcg.fits").astype(numpy.float32)
    e = ellipse.ellipse(dic,data,False, 1000.0)
    print("Read "+testingImg+" image")
    print(e.posx,e.posy)

    # Testing spirality profile:
    #testDistances(img,e,1000,numpy.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8],dtype=numpy.float32))

    galaxyIO.runSextractor(pathSp,testingImgSp,0)
    dic, data = galaxyIO.readSextractorOutput("0.cat")
    dic,data = galaxyIO.filterSextractorData(img, dic,data)
    bcg = gio.readFITSIMG("0_bcg.fits").astype(numpy.float32)
    e = ellipse.ellipse(dic,data,False, 1000.0)
    print("Read "+testingImg+" image")
    print(e.posx,e.posy)

    # Testing spirality profile:
    testDistances(imgSp,e,1000,)
    #print("Result: ", par.spirality3(imgSp,400,e))
    
    #saveButterworthProfile(2)
    #saveButterworthProfile(4)
    #saveButterworthProfile(8)
    #saveButterworthProfile(16)


    #smoothingTest(img)
    #rotationTest(img)

    # Sampling test
    #print("Testing sampling at dist = 10")
    #ps1,ps2 = gridAlg.getManySample(img,e,100,350)
    #print(ps1)
    #print("A",ps1[0])
    #scatterPlot(ps1[0],ps1[1],0, int(e.posx)-12, int(e.posy)+12,'b')
    #scatterPlot(ps2[0],ps2[1],0, int(e.posx)-12, int(e.posy)+12,'r')
    #plt.gca().add_patch(plt.Circle((e.posx,e.posy), radius=10, fill=False,linewidth=2.5))
    #plt.axes().arrow(e.posx,e.posy, 7.0, 6.0, head_width=0.5,head_length=0.5,fc='k', ec='k')
    #plt.text(123.0, 123.5,'R')
    #plt.show()

    #Concentration profile test:
    #print("Testing concentration profile...")
    #concentrationTotalTest(img,e,50,15)
    #print("Done")

    #Entropy test:
    #print("Testing Entropy Max-Min transition")
    #entropyTestMaxMin("tests/entropy/testeEntropiasExtremos.txt")


    #print("Testing K Folds (Entropy)")
    #kfolds = [i for i in range(2,1000,1)]
    #ellipticalEntropies = entropyBinTest(ellipticaWN,kfolds)
    #print("Starting entropy kbins")
    #spiralEntropies = entropyBinTest(spiralWN,kfolds)
    #normEntropyDiff = numpy.array([ellipticalEntropies[i]-spiralEntropies[i] for i in range(len(ellipticalEntropies))])
    #numpy.savetxt('tests/entropy/kBins.txt',[kfolds,normEntropyDiff])
    #print("Entropy kbins done!")
    
    #Sky box test:
    #print("Testing SkyBox (concentration)")
    #skyBoxTest(img)

    # Sampling Error test:
    #print("Testing sampling (concentration)")
    #samplingErrorTest("tests/concentration/erroAmostragem50.txt",density=50)
    #samplingErrorTest("tests/concentration/erroAmostragem100.txt",density=100)
    #samplingErrorTest("tests/concentration/erroAmostragem150.txt",density=150)
    
    plt.show()
