import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
#from scipy.ndimage import gaussian_filter
#import time

plt.close('all') #chiude tutte le figure aperte
#probabilmente devo avere Npxtip<280 circa (16gb ram con 25% circa occupato)
#278 funziona, 296 no
Npx=1500 #resolution of the map
pxlen=2 #nm, phyisical length of a pixel in the map
thres=10 #nm
s=0.5 #scaling factor to eventually avoid memory error (dilation)

rspikemin=s*50
rspikemax=s*90
hspikemin=s*150
hspikemax=s*200

rtip=s*300
#rtip=s*500

htip=hspikemax*1.02
spikedist=2*np.sqrt(rtip**2 - (rtip-hspikemax)**2) + 2
#MAP-------------------------------------------------
# z=mf.genFlat(Npx)
# #z=mf.genNormNoise(z,pxlen, 2, 2)
# z=mf.genHexSpikes(z,pxlen,hspikemin,hspikemax,spikedist,rspikemin,rspikemax,
#                   0.9,0.2,par='r', xmin=spikedist/2, xmax=len(z)*pxlen-spikedist/2,
#                   ymin=spikedist/2, ymax=len(z)*pxlen-spikedist/2)

# mf.plotfalsecol(z,pxlen)
# np.savetxt('map1.dat',z, header=str(pxlen)+' '+str(Npx)+' '+str(rtip)+' '+str(s)+' \n pxlen, Npx, Rtip, s')
# #Tip------------------------------------------------
# #occhio che h/a>>pxlen
# if hspikemax>rtip: print('Warning: possible spikes higher than semisphere tip')
# tip=mf.genSemisphTip(pxlen,htip,r=rtip)
# np.savetxt('tip1.dat',tip, header=str(pxlen)+' '+str(Npx)+' '+str(rtip)+' '+str(s)+' \n pxlen, Npx, Rtip, s')
# mf.plotfalsecol(tip,pxlen)
# # #IMG------------------------------------------------
# img = mph.grey_dilation(z, structure=-tip)
# img = mf.genNormNoise(img,pxlen, 4, 2)
# np.savetxt('img1.dat',img, header=str(pxlen)+' '+str(Npx)+' '+str(rtip)+' '+str(s)+' \n pxlen, Npx, Rtip, s')
# mf.plotfalsecol(img,pxlen)
# #---------------------------------------------------

img = np.loadtxt('img1.dat')
img = mph.filters.median_filter(img, 5)
img_obj_list, img_labeled, img_obj_ind = mf.identObj(img, thres, Npx_min=5)
mf.plotThres(img, img_labeled, pxlen)