import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
import part_num

plt.close('all') #chiude tutte le figure aperte

Npx=4000 #resolution of the map
pxlen=1 #nm, phyisical length of a pixel in the map
rmin=10 #nm
rmax=30 #nm
thres=0 #nm
Npart=100 #number of particles on the map
h = 80
r = 4

#MAP-------------------------------------------------

z=mf.genFlat(Npx)
#z=mf.genHexSpikes(z,pxlen,20,2,20)
#z=mf.genNormNoise(z,pxlen,100,50)
#z=mf.genFlat(Npx)
#z=gaussian_filter(z,20)
#z=mf.genSphere(z,pxlen,np.array([100,100,20,70]),np.array([10,30]))
#z=mf.genUnifSph(z,pxlen,10,5,10) #, xmin=22, xmax=58, ymin=62, ymax=78)
z=mf.genNormSph(z,pxlen,Npart,10,2) #, xmin=22, xmax=58, ymin=62, ymax=78)
#z=mf.genUnifIsolSph(z,pxlen,Npart,rmin,rmax) #, xmin=22, xmax=58, ymin=62, ymax=78)

#Tip------------------------------------------------
#occhio che h/a>>pxlen

tip=mf.genParabolicTip(pxlen,h,r=r) 
#tip=mf.genPyramidTip(pxlen,50,angle=np.pi/4)
#tip=mf.genConeTip(pxlen,50,angle=np.pi/2)
#tip=mf.genSemisphTip(pxlen,50,r=20)

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
mf.plotview(tip,pxlen,30,30)
#mf.plotview(tip,pxlen,90,0)
#mf.plotview(tip,pxlen,0,90)

#PARTNUM--------------------------------------------

part_num.reconstLogNorm(z, pxlen, thres, Npart)