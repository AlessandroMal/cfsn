import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from scipy.ndimage import gaussian_filter
#import time

plt.close('all') #chiude tutte le figure aperte

Npx=500 #resolution of the map
pxlen=1 #nm, phyisical length of a pixel in the map
thres=3 #nm

rspikemin=20 #nm
rspikemax=40 #nm
hspikemin=20
hspikemax=30

htip=120
ripmin=40
rtipmax=60

spikedist=100
#MAP-------------------------------------------------

z=mf.genFlat(Npx)
z=mf.genHexSpikes(z,pxlen,hspikemin,hspikemax,spikedist,rspikemin,rspikemax,0.8,0,
                  par='r', xmin=rtipmax, xmax=len(z)*pxlen-rtipmax,
                  ymin=rtipmax, ymax=len(z)*pxlen-rtipmax)
#DEVO
#OVALIZZARE LE SPIKES UNIFORMEMENTE
#FARE FILTRO
#Tip------------------------------------------------
#occhio che h/a>>pxlen
tip=mf.genRndSemisphTip(pxlen,120,40,60)

#IMG------------------------------------------------

img = mph.grey_dilation(z, structure=-tip)

#PLOT-----------------------------------------------

mf.plotfalsecol(z,pxlen)
#mf.plotview(z,pxlen,90,-90)
mf.plotfalsecol(img,pxlen)
#mf.plotfalsecol(tip,pxlen)

obj = mf.identObj(img,thres)
#for i in obj: mf.plotfalsecol(i,pxlen)
filtobj= par.revimg_filter(obj,pxlen,thres,0.5,1.5,0.2)
#for i in obj: print(par.capPar(i,pxlen,thres))
#for i in filtobj: mf.plotfalsecol(i,pxlen)

def V_h(h_obj, R): return np.pi/3*h_obj**2 *(3*R-h_obj)
def A_h(h_obj, R): return np.pi*h_obj *(2*R-h_obj)
def V_A(A_obj, R): return 2*np.pi*R**3 /3 *(1 - np.sqrt(1-A_obj/R**2 /np.pi) - A_obj/R**2 /np.pi/2) *(2 + np.sqrt(1-A_obj/R**2 /np.pi))

x_par=[]
y_par=[]
for i in filtobj:
    x_par.append(par.capPar(i,pxlen,thres)[0]) #h è 0, A è 2, V è 3
    y_par.append(par.capPar(i,pxlen,thres)[3])
#x_par=np.array(x_par)
#y_par=np.array(y_par)

R_opt, R_var= curve_fit(V_h, x_par, y_par)
print('Rtip=',*R_opt)
print('Rvar=',np.sqrt(*R_var))