from libc.math cimport log10,sqrt, atan2, pow, fabs
import numpy
import gridAlg
import scipy.ndimage as ndimage
import scipy.signal as signal
import math
from scipy.optimize import curve_fit
from scipy import fftpack
from random import shuffle

cimport numpy
cimport cython
cdef float pi = 3.14159265

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef float[:] ravel(float[:,:] mat):
    cdef:
          int w, h, countNotMasked, j, i,it
          float[:] line
    w, h = len(mat[0]), len(mat)
    countNotMasked = 0
    for i in range(w):
        for j in range(h):
            if(mat[j, i] != 0.0):
               countNotMasked = countNotMasked + 1
    line = numpy.array([0.0 for i in range(countNotMasked)], dtype=numpy.float32)
    it = 0
    for i in range(w):
        for j in range(h):
            if(mat[j, i] != 0.0):
                line[it] = mat[j, i]
                it = it + 1 
    return line
     

############################################
#     Funcoes de entropia:
# entrada:
#    mat - matriz do tipo numpy(list(list))
#    bins - numero de patamares
# retorna:
#    coeficiente encontrado,
#    matriz
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def entropyFunction(float[:,:] mat,int bins):
    cdef:
        float[:] freq, line
        double[:] binagem
        long[:] temp
        float somatorio,coef
        tuple x
        list entropies
        int w, h, i
    w, h = len(mat[0]), len(mat)
    line = ravel(mat)
    freq = numpy.array([0.0 for i in range(bins)], dtype=numpy.float32)
    temp, binagem = numpy.histogram(line,bins)
    somatorio = 0.0
    for i in range(bins):
        somatorio = somatorio + temp[i]
    for i in range(bins):
        freq[i] = float(temp[i])/float(somatorio)
    somatorio = 0.0
    for i in range(bins):
        if freq[i]>0.0:
            somatorio = somatorio - freq[i]*log10(freq[i])
    coef = somatorio/log10(bins)
    return coef


############################################
#     Translating and rotating
def rotateImage(float[:,:] img,float angle):
    imgR = ndimage.rotate(img, angle, reshape=False,mode='nearest')
    return imgR

############################################
#     Funcoes de assimetria:
# entrada:
#    mat - matriz do tipo numpy(list(list))
#    corrFunction - funcao de correlacao (stats.pearsonr, stats.spearmanr, ....)
# retorna:
#    coeficiente encontrado,
#    matrizRotacionada
#    pontos de correlacao (I,I_h)
#    mascara de pontos considerados (rotacionado e nao rotacionado)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def asymmetryFunction(float[:,:] mat,float[:,:] mask, corrFunct, ell):
    cdef:
         float minimo
         float[:,:]  matInverse, maskInverse,matWBCG
         int w, h, countNotMasked, it
    w, h = len(mat[0]), len(mat)
 
    maskInverse = rotateImage(mask, 90.0)
    matInverse = rotateImage(mat, 90.0)
    
    countNotMasked = 0
    for i in range(w):
        for j in range(h):
            if (mask[j,i] < 0.5) and (maskInverse[j,i] < 0.5)and (sqrt((i-w/2)**2+(j-h/2)**2) > 0.2*ell.fwhm): 
                countNotMasked = countNotMasked + 1
    v1 = []
    v2 = []
    it = 0
    for i in range(w):
        for j in range(h):
            if (mask[j,i] <= 0.5) and (maskInverse[j,i] <= 0.5)and (sqrt((i-w/2)**2+(j-h/2)**2) > 0.2*ell.fwhm): 
                v1.append(mat[j,i])
                v2.append(matInverse[j,i])
                it = it + 1

    mv1 = numpy.max(v1)
    mv2 = numpy.max(v2)
    for it in range(countNotMasked):
        v1[it] = v1[it]/mv1
        v2[it] = v2[it]/mv2
    

    coef = 1.0-corrFunct(v1, v2)[0]

    return coef, matInverse, (v1, v2)

#     Funcao Concentration :
############################################
# entrada:
#    mat - matriz do tipo numpy(list(list))
#    calibratedData - booleano indicando se os dados estao calibrados
#    p1 - porcentagem de luminosidade do numerador (0.8->c1,  0.9 -> c2,)
#    p2 - porcentagem de luminosidade do denominador (0.2 -> c1,  0.5 -> c2)
#    ell - elipse (com o centro definido)
#    nbins - variacao do raio avaliado
# retorna:
#    coeficiente encontrado,
#    matriz,
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def concentrationFunction(dists,concSeq, p1, p2, minCut):
    total, cutDist = getTotalLum(dists,concSeq,1.5,minCut)
    rn = gridAlg.findRadiusLuminosity(dists,concSeq,total, p1)
    rd = gridAlg.findRadiusLuminosity(dists,concSeq,total, p2)
    
    if(rd > 0.0):
        return log10(rn/rd),total
    else:
        raise Exception("Zero Division Error in concentration!")

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def concentrationFunction2(dists,concSeq, p1, p2, minCut):
    total = sum(concSeq)
    rn = gridAlg.findRadiusLuminosity(dists,concSeq,total, p1)
    rd = gridAlg.findRadiusLuminosity(dists,concSeq,total, p2)
    
    if(rd > 0.0):
        return log10(rn/rd),total
    else:
        raise Exception("Zero Division Error in concentration!")

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def getTotalLum(dists,conc, k,minD):

    # concentration central finite difference
    dconc = numpy.array([0.0 for i in range(len(conc)-2)])
    for i in range(1,len(conc)-1):
        dconc[i-1] = (conc[i+1] - conc [i-1]) / (dists[i+1]-dists [i-1])

    dconc = numpy.gradient(conc)
    # Finding distance that dconc < 1%
    for cutDist in range(2,len(dconc)):
        temp = conc[i] if  conc[i] != 0.0 else  1.0
        if(abs(dconc[cutDist]/conc[cutDist])<minD):
            break 

    if(len(dconc)-1 == cutDist):
        raise Exception("Not Convergent Concentration!")

    median = numpy.median(dconc[cutDist:len(dconc)])
    sigma = numpy.median(abs(dconc[cutDist:len(dconc)] - median))
    avg = 0.0
    n = 0.0
    for i in range(1,len(dconc)):
        if(abs(dconc[i])<k*sigma):
            n = n + 1.0
            avg += conc[i]
    if n > 10:
        return avg/n, cutDist
    else:
        raise Exception("Not Convergent Concentration!")


############################################
#     Funcao Smoothness (Suavizacao):
# entrada:
#    mat - matriz do tipo numpy(list(list))
#    kernel - matriz de convolucao
#    corrFunction - funcao de correlacao (stats.pearsonr, stats.spearmanr, ....)
# retorna:
#    coeficiente encontrado (valor e prob. da hipotese nula),
#    matriz,
#    matrizRotacionada
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def smoothnessFunction(float[:,:] mat, float[:,:]smMat, float[:,:] mask, corrFunct, ell):
    cdef:
        int w, h, countPts,it,i,j,countNotMasked
    w, h = len(mat[0]), len(mat)

    
    # counting the number of segmented pixels 
    countPts=0
    for i in range(w):
        for j in range(h):
            if(mask[j, i] == 0.0) and (sqrt((i-w/2)**2+(j-h/2)**2) > 0.2*ell.fwhm):
                countPts += 1

    
    if(countPts<6):
        raise Exception("Invalid number of smoothing pixels")
    it = 0 
    
    v1,v2= [0.0 for i in range(countPts)],[0.0 for i in range(countPts)]

    for i in range(w):
        for j in range(h):
            if (mask[j,i] == 0.0)and (sqrt((i-w/2)**2+(j-h/2)**2) > 0.2*ell.fwhm): 
                v1[it] = mat[j,i]
                v2[it] = smMat[j,i]
                it = it + 1
                
    #countNotMasked = len(v1)
    #mv1 = numpy.max(v1)
    #mv2 = numpy.max(v2)
    #for it in range(countNotMasked):
    #    v1[it] = v1[it]/mv1
    #    v2[it] = v2[it]/mv2

    coef = 1.0 - corrFunct(v1, v2)[0]

    return coef,(v1, v2)
############################################
#     Funcao Smoothness (Suavizacao):
# entrada:
#    mat - matriz do tipo numpy(list(list))
#    kernel - matriz de convolucao
#    corrFunction - funcao de correlacao (stats.pearsonr, stats.spearmanr, ....)
# retorna:
#    coeficiente encontrado (valor e prob. da hipotese nula),
#    matriz,
#    matrizRotacionada
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def smoothness1(float[:,:] mat, float[:,:]smMat, float[:,:] mask,ell):
    cdef:
        int w, h, countPts,it,i,j,countNotMasked
        float sm
    w, h = len(mat[0]), len(mat)

    
    # counting the number of segmented pixels 
    countPts=0
    for i in range(w):
        for j in range(h):
            if(mask[j, i] == 0.0) and (sqrt((i-w/2)**2+(j-h/2)**2) > 0.2*ell.fwhm):
                countPts += 1

    
    if(countPts<6):
        raise Exception("Invalid number of smoothing pixels")
    sm = 0.0
    sm2 = 0.0
    for i in range(w):
        for j in range(h):
            if (mask[j,i] == 0.0)and (sqrt((i-w/2)**2+(j-h/2)**2) > 0.2*ell.fwhm): 
                sm = abs(mat[j,i]-smMat[j,i])
                sm2 += abs(mat[j,i])
    return sm/sm2


############################################
#     Funcao Spirality (Espiralidade):
# entrada:
#    mat - matriz do tipo numpy(list(list))
#    nPhase - numero de angulos na matriz  em coordenadas polares
#    nRadius - numero de distancias na matriz  em coordenadas polares
#    ellipse - elipse que contem galaxia
# retorna:
#    coeficiente encontrado ,
#    matriz,
#    matriz em Coordenadas Polares
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def spiralityFunction(float[:,:] mat,float[:,:]  mask, int nPhases,int  nRadius, ellipse):
    cdef:
        float[:,:] transformed, transformedMask,x, y
        float phase, mod, spirality, sumPhase, nphase
        int yIt, xIt
    transformed = gridAlg.transformPolarCoord(mat, ellipse, nPhases, nRadius)
    transformedMask = gridAlg.transformPolarCoord(mask, ellipse, nPhases, nRadius)
    y, x = gridAlg.gradient(transformed)

    ## y -> phase
    ## x -> module
    sumPhase = float(0.0)
    nphase = 1
    for yIt in range(len(y)):
        for xIt in range(len(y[yIt])):
            if(transformedMask[yIt, xIt] >= 0.5):
                continue
            mod = sqrt(pow(y[yIt, xIt],2.0)+pow(x[yIt, xIt],2.0))
            phase = atan2(y[yIt, xIt],x[yIt, xIt])
            phase = phase if phase > 0.0 else phase + 2.0*pi
            y[yIt, xIt] = phase
            x[yIt, xIt] = mod
            sumPhase = sumPhase + phase
            nphase = nphase + 1
    if(nphase > 0):
        mean = sumPhase / float(nphase)
    else:
        raise Exception("No points found in transformed polar coordinates masked!")
    spirality = 0.0
    for yIt in range(len(y)):
        for xIt in range(len(y[yIt])):
            if(transformedMask[yIt, xIt] >= 0.5):
                continue
            spirality += pow(y[yIt, xIt] - mean,2.0)
    if(nphase > 1):
        spirality = sqrt(spirality / float(nphase-1))
    else:
        raise Exception("No points found in transformed polar coordinates masked!")
    return (spirality, transformed, x, y)

############################################
#     Funcao Spirality (Espiralidade):
# entrada:
#    mat - matriz do tipo numpy(list(list))
#    nPhase - numero de angulos na matriz  em coordenadas polares
#    nRadius - numero de distancias na matriz  em coordenadas polares
#    ellipse - elipse que contem galaxia
# retorna:
#    coeficiente encontrado ,
#    matriz,
#    matriz em Coordenadas Polares
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def spiralityEllipsoidalFunction(float[:,:] mat,float[:,:]  mask,int nPhases,int nRadius, ellipse):
    cdef:
        float[:,:] transformed, transformedMask, x, y
        float phase, mod, spirality, sumPhase, nphase
        int yIt, xIt
    transformed = gridAlg.transformEllipseCoord(mat, ellipse, nPhases, nRadius)
    transformedMask = gridAlg.transformEllipseCoord(mask, ellipse, nPhases, nRadius)
    
    y, x = gridAlg.gradient(transformed)

    ## y -> phase
    ## x -> module
    sumPhase = float(0.0)
    nphase = 1
    for yIt in range(len(y)):
        for xIt in range(len(y[yIt])):
            if(transformedMask[yIt, xIt] >= 0.5):
                continue
            mod = sqrt(pow(y[yIt, xIt],2.0)+pow(x[yIt, xIt],2.0))
            phase = atan2(y[yIt, xIt],x[yIt, xIt])
            phase = phase if phase > 0.0 else phase + 2.0*pi
            y[yIt, xIt] = phase
            x[yIt, xIt] = mod
            sumPhase = sumPhase + phase
            nphase = nphase + 1
    if(nphase > 0):
        mean = sumPhase / float(nphase)
    else:
        raise Exception("No points found in transformed spiral coordinates masked!")
    spirality = 0.0
    for yIt in range(len(y)):
        for xIt in range(len(y[yIt])):
            if(transformedMask[yIt, xIt] >= 0.5):
                continue
            spirality += pow(y[yIt, xIt] - mean,2.0)
    if(nphase > 1):
        spirality = sqrt(spirality / float(nphase-1))
    else:
        raise Exception("No points found in transformed spiral coordinates masked!")
    return (spirality, transformed, x, y)
    


def spirality3(float[:,:] mat,float[:,:] bcg, ellipse):
    tseq = []
    nsamplePerDist=1000
    distances=numpy.array([0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75],dtype=numpy.float32)
    mat2 = numpy.array([[ mat[i][j]-bcg[i][j] for i in range(len(mat[i]))] for j in range(len(mat))], dtype=numpy.float32)
    for d in distances:
        seq, trash = gridAlg.sampleSameDistance(mat2,nsamplePerDist,d,ellipse)
        avg = numpy.average(seq)
        fft = fftpack.fft(seq/avg)
        tseq=numpy.concatenate([tseq,numpy.real(fft[1:(len(fft)-1)])])
        tseq=numpy.concatenate([tseq,numpy.imag(fft[1:(len(fft)-1)])])
    return max(numpy.abs(tseq))
    

############################################
#     Funcao Hamming kernel:
# entrada:
#    dim - dimensao da matriz
# retorna:
#    matriz resultante
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def hammingKernel(int dim):
    line = signal.hamming(float(dim))
    mat = numpy.sqrt(numpy.outer(line,line))
    mat = mat /mat.sum()
    mat = (mat).astype('float32')
    return mat

############################################
#     Funcao box car kernel:
# entrada:
#    dim - dimensao da matriz
# retorna:
#    matriz resultante
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def boxCarKernel(int dim):
    mat = numpy.array([[1.0 for i in range(dim)] for j in range(dim)])
    mat = mat / mat.sum()
    return mat



