import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
#from scipy.ndimage import gaussian_filter
#import time

# import random
# random.seed(1)
# np.random.seed(1)

plt.close('all') #chiude tutte le figure aperte

Npx=2000 #resolution of the map
pxlen=1 #nm, phyisical length of a pixel in the map
rmin=20 #nm
rmax=20 #nm
thres=15 #nm
Npart=250 #number of particles on the map
h = 80
ar = 2

#MAP-------------------------------------------------

z=mf.genFlat(Npx)
#z=mf.genHexSpikes(z,pxlen,20,2,20)
#z=mf.genNormNoise(z,pxlen,100,50)
#z=mf.genFlat(Npx)
#z=gaussian_filter(z,20)
#z=mf.genSphere(z,pxlen,np.array([100,100,20,70]),np.array([10,30]))
#z=mf.genUnifSph(z,pxlen,10,5,10) #, xmin=22, xmax=58, ymin=62, ymax=78)
#z=mf.genNormSph(z,pxlen,Npart,20,0) #, xmin=22, xmax=58, ymin=62, ymax=78)
z=mf.genUnifIsolSph(z,pxlen,Npart,rmin,rmax) #, xmin=22, xmax=58, ymin=62, ymax=78)

#Tip------------------------------------------------
#occhio che h/a>>pxlen

tip=mf.genParabolicTip(pxlen,h,ar) 
#tip=mf.genPyramidTip(pxlen,50,2)
#----------la dilation non funziona bene con piramide ??
#tip=mf.genSemisphTip(pxlen,40,1)

#IMG------------------------------------------------

img = mph.grey_dilation(z, structure=-tip)

#PLOT-----------------------------------------------

mf.plotview(z,pxlen,30,30)
#mf.plotview(z,pxlen,90,0)
#mf.plotview(z,pxlen,0,90)
mf.plotfalsecol(z,pxlen)
mf.plotview(img,pxlen,30,30)
#mf.plotview(img,pxlen,90,0)
#mf.plotview(img,pxlen,0,90)
mf.plotfalsecol(img,pxlen)
#mf.plotview(tip,pxlen,30,30)
#mf.plotview(tip,pxlen,90,0)
#mf.plotview(tip,pxlen,0,90)

#TIME---------------------------------------------
'''
v=[]
a=100
b=1000
step=100
b+=1
for i in range(a,b,step):
    z=mf.genFlat(i)
    start=time.perf_counter()
    par.specArea(z,pxlen)
    v.append(time.perf_counter()-start)
    
plt.figure(figsize=(20,20))
plt.plot(np.arange(a,b,step), v, marker='.')
plt.plot(np.arange(a,b,step),6.5/9/10000*np.arange(a,b,step)**2)
plt.grid()
'''  
#PARAMS---------------------------------------------
    
#print(8*2/3*3.14*20**3+8*3.14*20**3,"  ",par.V(z,pxlen),"  ",par.V(img,pxlen))
#zParams = par.calcParams(z,pxlen,thres)
#print('\nSURFACE:\n', zParams)
#print(par.specArea(z,pxlen))
#print(time.perf_counter()-start)
#imgParams = par.calcParams(img,pxlen,thres)
#print('\nIMAGE:\n', imgParams)

#par.paramvsTip(z, pxlen, thres, mf.genParabolicTip, 80, 0.5, 10, 5)
#par.paramvsNpart(Npx, pxlen, rmin, rmax, 1, 20,
                # mf.genParabolicTip, 80, 0.5, 5, 4,
                # 1, lambda surf: par.coverage(surf, thres), r'$coverage$')
# par.paramvsRpart(Npx, pxlen, 1, 64, 10, 5,
#                  mf.genParabolicTip, 80, 0.5, 5, 4,
#                  1, np.mean, r'$mean [nm]$')

par.singlePartAnalysis(z, img, pxlen, thres, h, ar, Npart, Npx, 
                        param_string_list = ['max [nm]', 'mean [nm]', 'std[nm]', 'skewness', 'V[nm³]', 'specArea', 'coverage'],
                        binwidth_list=[0.1, 0.25, 0.25, 0.2, 5000, 0.25, 0.025],
                        bin_min_list=[None,None,None,None,None,None,None],
                        bin_max_list=[None,None,None,None,None,None,None])
plt.show()
#[0.1, 2.5, 1.5, 0.2, 5000, 0.25, 0.025]
#[0.1, 0.25, 0.25, 0.2, 5000, 0.25, 0.025]