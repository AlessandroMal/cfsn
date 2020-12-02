import mapsfunction as mf
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt

plt.close('all') #chiude tutte le figure aperte

Npx=300
pxlen=1 #nm

#MAP-------------------------------------------------

z=mf.genFlat(Npx)
#z=mf.genNormNoise(z,pxlen,2,50)
#z=mf.genSphere(z,pxlen,np.array([100,100,20,70]),np.array([10,30]))
#z=mf.genUnifSph(z,pxlen,10,5,10) #, xmin=22, xmax=58, ymin=62, ymax=78)
#z=mf.genNormSph(z,pxlen,8,15,5) #, xmin=22, xmax=58, ymin=62, ymax=78)

#print(mf.specArea(z,pxlen))

z=mf.genUnifIsolSph(z,pxlen,8,20,20) #, xmin=22, xmax=58, ymin=62, ymax=78)
#OCCHIO AI DATI IN INPUT PERCHE FUNZIONA A RIGETTO

#TIP------------------------------------------------

tip=mf.genParabolicTip(pxlen,80,4) #occhio che h/a>>pxlen
#tip=mf.genPyramidTip(pxlen,50,2) #occhio che h/a>>pxlen
#----------la dilation non funziona bene con piramide ??
#tip=mf.genSemisphTip(pxlen,40,1)

#IMG------------------------------------------------

img = mph.grey_dilation(z, structure=-tip)
#---------------------------------------------------

#print(8*2/3*3.14*20**3+8*3.14*20**3,"  ",mf.V(z,pxlen),"  ",mf.V(img,pxlen))

#PLOT-----------------------------------------------
mf.plotview(z,pxlen,30,30)
#mf.plotview(z,pxlen,90,0)
#mf.plotview(z,pxlen,0,90)
mf.plotfalsecol(z,pxlen)
#mf.plotview(img,pxlen,90,0)
mf.plotview(img,pxlen,30,30)
mf.plotview(img,pxlen,0,90)
mf.plotfalsecol(img,pxlen)
#mf.plotview(tip,pxlen,30,30)
#mf.plotview(tip,pxlen,90,0)
#mf.plotview(tip,pxlen,0,90)

#PARAMS---------------------------------------------

zParams = mf.calcParams(z,pxlen)
print('\nSURFACE:\n', zParams)
imgParams = mf.calcParams(img,pxlen)
print('\nIMAGE:\n', imgParams)

#mf.paramTipDepend(z, pxlen, 'genParabolicTip', 80, 1, 40, 20)

#obj = mf.identObj(z,0)
#for i in obj: mf.plotfalsecol(i,pxlen)

