from math import sin,cos,tan,pi,floor,log10,sqrt,pow,radians, fabs
import numpy as np
import pyfits
import copy
from astropy import wcs
import sys
import os
from subprocess import call
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
import ConfigParser

def replaceAfterLastDot(fileName, str2Rep):
    splitted = fileName.split('.')
    splitted[len(splitted)-1] = str2Rep
    return ".".join(splitted)

#inputs:
#    python maskmaker_wcut.py filename parallelIndex ra dec stampsize
def main():
    item = sys.argv[1]
    pIndex = sys.argv[2]
    ra = float(sys.argv[3])
    dec = float(sys.argv[4])
    stampSize = float(sys.argv[5])    

    infilename = item#item[0:len(item) - 1]
    ptsInside2Rp = []
    print("File: " + item)
    header = pyfits.getheader(infilename, 0)
    data = pyfits.getdata(infilename, 0)
    ylen, xlen = data.shape[0], data.shape[1]
    print("Tamanho:", ylen,xlen)
    CRPIX1 = float(header['CRPIX1'])
    CRPIX2 = float(header['CRPIX2'])
    CRVAL1 = float(header['CRVAL1'])
    CRVAL2 = float(header['CRVAL2'])
    CD1_1 = float(header['CD1_1'])
    CD1_2 = float(header['CD1_2'])
    CD2_1 = float(header['CD2_1'])
    CD2_2 = float(header['CD2_2'])
    if header['CTYPE1'] == 'DEC--TAN':
        CD1_1 = float(header['CD2_1'])
        CD2_1 = float(header['CD1_1'])
        CD1_2 = float(header['CD2_2'])
        CD2_2 = float(header['CD1_2'])
        CRVAL1 = float(header['CRVAL2'])
        CRVAL2 = float(header['CRVAL1'])
    w = wcs.WCS(infilename, relax=True)
    w.wcs.crpix = [CRPIX1, CRPIX2]
    w.wcs.crval = [CRVAL1, CRVAL2]
    w.wcs.cd = [[CD1_1, CD2_1], [CD1_2, CD2_2]]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    sc = SkyCoord(ra=ra, dec=dec, unit='deg')
    px, py = skycoord_to_pixel(sc,w,0)
    objcoord = [px, py]

    with open('base_default.sex', 'r') as inpsx:
        with open('default'+pIndex+'.sex', 'w') as outsx:
            for item2 in inpsx.readlines():
                ll2 = item2.split()
                if len(ll2) > 0:
                    if ll2[0] == 'CATALOG_NAME':
                        ll2[1] = 'out_sex_large'+pIndex+'.cat'
                        lstrin = ' '
                        for k in range(0, len(ll2)):
                            lstrin += ll2[k] + ' '
                        outsx.write('%s\n' % lstrin[1:len(lstrin)])
                    elif ll2[0] == 'BACK_SIZE':
                        ll2[1] = '128'
                        lstrin = ' '
                        for k in range(0, len(ll2)):
                            lstrin += ll2[k] + ' '
                        outsx.write('%s\n' % lstrin[1:len(lstrin)])
                    elif ll2[0] == 'CHECKIMAGE_NAME':
                        ll2[1] = 'check1_large.fits, check2_large.fits, check3_large.fits'
                        lstrin = ' '
                        for k in range(0, len(ll2)):
                            lstrin += ll2[k] + ' '
                        outsx.write('%s\n' % lstrin[1:len(lstrin)])
                    else:
                        outsx.write('%s' % item2)
                else:
                    outsx.write('%s' % item2)

        configFile = ConfigParser.ConfigParser()
        configFile.read('../cfg/paths.ini')
        sexPath = configFile.get("Path","Sextractor")
        call(sexPath+' '+infilename+' -c default'+pIndex+'.sex',shell=True)

        # finding point:
        with open('out_sex_large'+pIndex+'.cat') as inflarge:
            dist=1000.
            for item2 in inflarge.readlines():
                ll2=item2.split()
                if ll2[0]!='#':
                        X=float(ll2[7])
                        Y=float(ll2[8])
                        XI=(CD1_1*(X-CRPIX1)+CD1_2*(Y-CRPIX2))*pi/180.
                        ETA=(CD2_1*(X-CRPIX1)+CD2_2*(Y-CRPIX2))*pi/180.
                        p=(XI**2+ETA**2)**0.5
                        c=np.arctan(p)
                        RA=np.arctan(XI*sin(c)/(p*cos(CRVAL2*pi/180.)*cos(c)-ETA*sin(CRVAL2*pi/180.)*sin(c)))*180./pi+CRVAL1
                        DEC=np.arcsin(cos(c)*sin(CRVAL2*pi/180.)+ETA*sin(c)*cos(CRVAL2*pi/180.)/p)*180./pi
                        if ((X-objcoord[0])**2+(Y-objcoord[1])**2)**0.5<dist:
                            dist=((X-objcoord[0])**2+(Y-objcoord[1])**2)**0.5
                            ARA=RA
                            ADEC=DEC
                            AX=int(X)
                            AY=int(Y)
                            APETROSIAN = float(ll2[4])*float(ll2[9])
                            thetabcg=float(ll2[11])
                            siz=float(ll2[9])*float(ll2[4])
                            razax=float(ll2[10])/float(ll2[9])
                            mag=float(ll2[2])
                            countbcg=float(ll2[1])

        #cutting image:
	
        #objcoord[0] = objcoord[0]-np.max([1,AX-siz])
        #objcoord[1] = objcoord[1]-np.max([1,AY-siz])
        #data = data[np.max([1,AY-siz]):np.min([data.shape[0],AY+siz]),np.max([1,AX-siz]):np.min([data.shape[1],AX+siz])]
        oldSiz = siz
        sizy = np.min(np.array([fabs(AY-siz),stampSize*fabs(siz), fabs(data.shape[0]-AY)]))
        sizx = np.min(np.array([fabs(AX-siz),stampSize*fabs(siz), fabs(data.shape[1]-AX)]))  
        siz = int(np.min(np.array([sizy,sizx])))
        data = data[AY-siz:AY+siz,AX-siz:AX+siz]  
        os.remove('out_sex_large'+pIndex+'.cat')
        os.remove('default'+pIndex+'.sex')
        if(siz != int(stampSize*fabs(oldSiz)) ):
             error = 2
        else:
             error = 0
        maskFileName = '../cutted/'+str(pIndex)+"_mask.fit"
        infilename = '../cutted/'+str(pIndex)+".fit"
        with open("../"+pIndex+"_log.txt","w") as l:
            l.write(str(error))
        pyfits.writeto('../cutted/'+str(pIndex)+".fit",data,clobber=True)       
