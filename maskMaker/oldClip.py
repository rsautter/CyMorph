from math import sin,cos,tan,pi,floor,log10,sqrt,pow,radians, fabs
import numpy as np
import pyfits
import copy
from astropy import wcs
import sys
from subprocess import call
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord

def replaceAfterLastDot(fileName, str2Rep):
    splitted = fileName.split('.')
    splitted[len(splitted)-1] = str2Rep
    return ".".join(splitted)

#inputs:
#    python maskmaker_wcut.py filename parallelIndex ra dec
#    python maskmaker_wcut.py filename parallelIndex -cart x y
def main():
    item = sys.argv[1]
    pIndex = sys.argv[2]
    ra = float(sys.argv[3])
    dec = float(sys.argv[4])
    
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
        with open('default.sex', 'w') as outsx:
            for item2 in inpsx.readlines():
                ll2 = item2.split()
                if len(ll2) > 0:
                    if ll2[0] == 'CATALOG_NAME':
                        ll2[1] = 'out_sex_large.cat'
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


        call('./../sex '+infilename+' -c default.sex',shell=True)

        # finding point:
        with open('out_sex_large.cat') as inflarge:
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
                            AX=X
                            AY=Y
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
        sizy = np.min(np.array([fabs(AY-siz),2.*fabs(siz), fabs(data.shape[0]-AY)])) 
        sizx = np.min(np.array([fabs(AX-siz),2.*fabs(siz), fabs(data.shape[1]-AX)]))  
        siz = np.min(np.array([sizy,sizx]))
        objcoord[0] = siz
        objcoord[1] = siz
        data = data[AY-siz:AY+siz,AX-siz:AX+siz]      

        maskFileName = '../cutted/'+str(pIndex)+"_mask.fit"
        infilename = '../cutted/'+str(pIndex)+".fit"
        #pyfits.writeto('../cutted/'+str(pIndex)+".fit",data,clobber=True)
        #pyfits.writeto('../cutted/t'+str(pIndex)+".fit",data,clobber=True)
        #call('./sex '+infilename+' -c default.sex',shell=True)

        # rerun with stamp
        with open('out_sex_large.cat') as inflarge:
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
                            AX=X
                            AY=Y
                            APETROSIAN = float(ll2[4])*float(ll2[9])
                            thetabcg=float(ll2[11])
                            siz=float(ll2[9])*float(ll2[4])
                            razax=float(ll2[10])/float(ll2[9])
                            mag=float(ll2[2])
                            countbcg=float(ll2[1])


        print("Point location (x,y): ("+str(AX)+", "+str(AY)+")")
        # start mask making
        mask = copy.deepcopy(data)
        mask[:,:]=0

        with open('out_sex_large.cat') as inflarge:
            for item2 in inflarge.readlines():
                ll2=item2.split()
                if ll2[0]!='#':
                    X = float(ll2[7])
                    Y = float(ll2[8])
                    PETROSIAN = float(ll2[4])*float(ll2[9])
                    if AX==X and AY==Y:
                        theta=float(ll2[11])*pi/180.
                        thresh=float(ll2[6])
                        if float(ll2[4])>float(ll2[3]) and float(ll2[4])>0:
                            aa=float(ll2[9])*float(ll2[4])
                            bb=float(ll2[10])*float(ll2[4])
                        elif float(ll2[3])>0:
                            aa=float(ll2[9])*float(ll2[3])*1.5
                            bb=float(ll2[10])*float(ll2[3])*1.5
                        else:
                            aa=2.0
                            bb=2.0
                        for ii in range(0,data.shape[1]):
                            for jj in range(0,data.shape[0]):
                                distp=((jj-Y)**2+(ii-X)**2)**0.5
                                if (ii-X)!=0:
                                    phi=np.arctan((jj-Y)/(ii-X))
                                else:
                                    phi=pi/2.
                                radius=aa*bb/((bb*cos(phi-theta))**2+(aa*sin(phi-theta))**2)**0.5
#                                radius2=2.*(thresh/(pi*aa*bb))**0.5*aa
#                                radius2=2.*(thresh/(pi*aa*bb))**0.5*radius
                                radius2=radius
                                if radius2<distp:
                                    mask[jj,ii]=1
                                else:
                                    ptsInside2Rp.append((jj,ii))
                        break
        print("N inside ellipse points: ",len(ptsInside2Rp))

        with open('out_sex_large.cat') as inflarge:
            for item2 in inflarge.readlines():
                ll2=item2.split()
                if ll2[0]!='#':
                    X = float(ll2[7])
                    Y = float(ll2[8])
                    PETROSIAN = float(ll2[4])*float(ll2[9])
                    if AX!=X and AY!=Y and sqrt(pow(AX-X,2.0)+pow(AY-Y,2.0)) < 4.0*(APETROSIAN + PETROSIAN):
                        thresh=float(ll2[6])
                        theta=float(ll2[11])*pi/180.
                        if float(ll2[4])>float(ll2[3]) and float(ll2[4])>0:
                            aa=1.*float(ll2[9])*float(ll2[4])
                            bb=1.*float(ll2[10])*float(ll2[4])
                        elif float(ll2[3])>0:
                            aa=1.*float(ll2[9])*float(ll2[4])*2.
                            bb=1.*float(ll2[10])*float(ll2[4])*2.
                        else:
                            aa=1.0
                            bb=1.0
                        for point in ptsInside2Rp:
                                jj=point[0]
                                ii=point[1]
                                distp=((jj-Y)**2+(ii-X)**2)**0.5
                                if (ii-X)!=0:
                                    phi=np.arctan((jj-Y)/(ii-X))
                                else:
                                    phi=pi/2.
                                if ((bb*cos(phi-theta))**2+(aa*sin(phi-theta))**2)**0.5 != 0.0:
                                    radius=aa*bb/((bb*cos(phi-theta))**2+(aa*sin(phi-theta))**2)**0.5
                                radius2=1.5*(thresh/(pi*aa*bb))**0.5*radius if (pi*aa*bb) != 0.0 else 1.5*(thresh)**0.5*radius
                                if radius>=distp and radius2>=distp:
                                    mask[jj,ii]=1

        # creating small-grained sky source catalogue
        with open('base_default.sex','r') as inpsx:
            with open('default.sex','w') as outsx:
                for item2 in inpsx.readlines():
                    ll2=item2.split()
                    if len(ll2)>0:
                        if ll2[0]=='CATALOG_NAME':
                            ll2[1]='out_sex_small.cat'
                            lstrin=' '
                            for k in range(0,len(ll2)):
                                lstrin+=ll2[k]+' '
                            outsx.write('%s\n' % lstrin[1:len(lstrin)])
                        elif ll2[0]=='BACK_SIZE':
                            ll2[1]=str(siz/5.)
                            lstrin=' '
                            for k in range(0,len(ll2)):
                                lstrin+=ll2[k]+' '
                            outsx.write('%s\n' % lstrin[1:len(lstrin)])
                        elif ll2[0]=='CHECKIMAGE_NAME':
                            ll2[1]='check1_small.fits, check2_small.fits, check3_small.fits'
                            lstrin=' '
                            for k in range(0,len(ll2)):
                                lstrin+=ll2[k]+' '
                            outsx.write('%s\n' % lstrin[1:len(lstrin)])
                        else:
                            outsx.write('%s' % item2)
                    else:
                        outsx.write('%s' % item2)


        call('./../sex '+infilename+' -c  default.sex',shell=True)


        # read and find bcg once again
        with open('out_sex_small.cat') as inflarge:
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
                        AX=X
                        AY=Y

        # add new sources to mask
        with open('out_sex_small.cat') as inflarge:
            for item2 in inflarge.readlines():
                ll2=item2.split()
                if ll2[0]!='#':
                    X=float(ll2[7])
                    Y=float(ll2[8])
                    PETROSIAN = float(ll2[4])*float(ll2[9])
                    if AX!=X and AY!=Y and sqrt(pow(AX-X,2.0)+pow(AY-Y,2.0)) < 4.0*(APETROSIAN + PETROSIAN):
                        theta=float(ll2[11])*pi/180.
                        thresh=float(ll2[6])
                        if float(ll2[4])>float(ll2[3]) and float(ll2[4])>0:
                            aa=1*float(ll2[9])*float(ll2[4])
                            bb=1*float(ll2[10])*float(ll2[4])
                        elif float(ll2[3])>0:
                            aa=1*float(ll2[9])*float(ll2[3])*2.
                            bb=1*float(ll2[10])*float(ll2[3])*2.
                        else:
                            aa=1.0
                            bb=1.0
                        for point in ptsInside2Rp:
                            jj=point[0]
                            ii=point[1]
                            distp=((jj-Y)**2+(ii-X)**2)**0.5
                            if (ii-X)!=0:
                                phi=np.arctan((jj-Y)/(ii-X))
                            else:
                                phi=pi/2.
                            radius=aa*bb/((bb*cos(phi-theta))**2+(aa*sin(phi-theta))**2)**0.5
                            radius2=1.5*(thresh/(pi*aa*bb))**0.5*radius if (pi*aa*bb) != 0.0 else 1.5*(thresh)**0.5*radius
                            if radius>=distp and radius2>=distp:
                                mask[jj,ii]=1
        # write mask
        print("Writting")
        pyfits.writeto('../cutted/'+str(pIndex)+".fit",data,clobber=True)
        #pyfits.writeto(str(pIndex)+"_m.fit",mask,clobber=True)
