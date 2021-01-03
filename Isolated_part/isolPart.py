import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
#from scipy.ndimage import gaussian_filter
#import time

plt.close('all') #chiude tutte le figure aperte

Npx=300 #resolution of the map
pxlen=1 #nm, phyisical length of a pixel in the map
rmin=20 #nm
rmax=20 #nm
thres=3 #nm
Npart=40 #number of particles on the map

#MAP-------------------------------------------------

#z=mf.genFlat(Npx)
#z=mf.genHexSpikes(z,pxlen,20,2,20)
#z=mf.genNormNoise(z,pxlen,100,50)
#z=mf.genFlat(Npx)
#z=gaussian_filter(z,20) #non funziona?
#z=mf.genSphere(z,pxlen,np.array([100,100,20,70]),np.array([10,30]))
#z=mf.genUnifSph(z,pxlen,10,5,10) #, xmin=22, xmax=58, ymin=62, ymax=78)
#z=mf.genNormSph(z,pxlen,8,15,5) #, xmin=22, xmax=58, ymin=62, ymax=78)
#z=mf.genHexSpikes(z,pxlen,40,1,40)
#z=mf.genUnifIsolSph(z,pxlen,Npart,rmin,rmax) #, xmin=22, xmax=58, ymin=62, ymax=78)
#Tip------------------------------------------------
#occhio che h/a>>pxlen

#tip=mf.genParabolicTip(pxlen,50,r=8) 
#tip=mf.genPyramidTip(pxlen,50,angle=np.pi/4)
#tip=mf.genConeTip(pxlen,50,angle=np.pi/2)
#tip=mf.genSemisphTip(pxlen,50,r=20)

#IMG------------------------------------------------

#img = mph.grey_dilation(z, structure=-tip)

#PLOT-----------------------------------------------

#mf.plotview(z,pxlen,30,30)
#mf.plotview(z,pxlen,90,-90)
#mf.plotview(z,pxlen,0,90)
#mf.plotfalsecol(z,pxlen)
#mf.plotview(img,pxlen,30,30)
#mf.plotview(img,pxlen,90,-90)
#mf.plotview(img,pxlen,0,90)
#mf.plotfalsecol(img,pxlen)
#mf.plotview(tip,pxlen,30,30)
#mf.plotview(tip,pxlen,90,-90)
#mf.plotview(tip,pxlen,0,90)
#mf.plotfalsecol(tip,pxlen)

#PARAMS---------------------------------------------

#print('hmax ',par.h_max(z,10))
#print(mf.Ncp(pxlen,Npx,rmax))
#zParams = par.calcParams(z,pxlen,thres)
#print('\nSURFACE:\n', zParams)
#print(par.specArea(z,pxlen))
#print(time.perf_counter()-start)
#imgParams = par.calcParams(img,pxlen,thres)
#print('\nIMAGE:\n', imgParams)

'''
z=mf.genFlat(Npx)
z=mf.genUnifIsolSph(z,pxlen,Npart,rmin,rmax)
x,y=par.paramvsTip(z, pxlen, thres, mf.genParabolicTip, 50, 2, 30, 10)
x[0]=x[0]/rmax
x[1]=x[1]+r'$/ \langle R_{part} \rangle$'
par.plotParams(x,y)
'''

#21 5
#par.paramvsNpart(Npx, pxlen, rmin, rmax, 0, 40, 3,
#                 mf.genParabolicTip, 85, 'r', 8, 40, 3,
#                 1, lambda surf: par.coverage(surf, thres), r'$rel. coverage$',0.04, True)

par.paramvsNpart(Npx, pxlen, rmin, rmax, 0, 4,
                 mf.genParabolicTip, 85, 8, 40, 3,
                 1, lambda surf: par.coverage(surf, thres), r'$rel. coverage$')


'''
par.paramvsRpart(Npx, pxlen, 10, 40, 3, 20,
                 mf.genParabolicTip, 50, 8, 40, 5,
                 1, lambda surf: par.coverage(surf, thres), r'$coverage$')
'''
#obj = mf.identObj(img,thres)
#for i in obj: mf.plotfalsecol(i,pxlen)
