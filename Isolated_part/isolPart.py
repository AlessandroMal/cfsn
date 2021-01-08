import mapsfunction as mf
import parameters as par
#import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
#from scipy.ndimage import gaussian_filter
#import time

plt.close('all') #chiude tutte le figure aperte

Npx=400 #resolution of the map
pxlen=0.75 #nm, phyisical length of a pixel in the map
rmin=20 #nm
rmax=20 #nm
thres=3 #nm
Npart=40 #number of particles on the map

par.paramvsNpart(Npx, pxlen, rmin, rmax, 1, 40, 14,
                 mf.genPyramidTip, 65, 'angle', np.pi/8, np.pi/2, 1,
                 1, lambda surf: par.V(surf, pxlen), 'VvsN.dat')

par.paramvsRpart(Npx, pxlen, 15, 10, 30, 11,
                 mf.genPyramidTip, 65, 'angle', np.pi/8, np.pi/2, 1,
                 1, lambda surf: par.V(surf, pxlen), 'VvsR.dat')

par.paramvsNpart(Npx, pxlen, rmin, rmax, 1, 40, 14,
                 mf.genPyramidTip, 65, 'angle', np.pi/8, np.pi/2, 1,
                 1, lambda surf: par.h_av(surf, pxlen), 'hvsN.dat')

par.paramvsRpart(Npx, pxlen, 15, 10, 30, 11,
                 mf.genPyramidTip, 65, 'angle', np.pi/8, np.pi/2, 1,
                 1, lambda surf: par.h_av(surf, pxlen), 'hvsR.dat')

par.paramvsNpart(Npx, pxlen, rmin, rmax, 1, 40, 14,
                 mf.genPyramidTip, 65, 'angle', np.pi/8, np.pi/2, 1,
                 1, lambda surf: par.h_std(surf, pxlen), 'stdvsN.dat')

par.paramvsRpart(Npx, pxlen, 15, 10, 30, 11,
                 mf.genPyramidTip, 65, 'angle', np.pi/8, np.pi/2, 1,
                 1, lambda surf: par.h_std(surf, pxlen), 'stdvsR.dat')

par.paramvsNpart(Npx, pxlen, rmin, rmax, 1, 40, 14,
                 mf.genPyramidTip, 65, 'angle', np.pi/8, np.pi/2, 1,
                 1, lambda surf: par.coverage(surf, thres), 'covvsN.dat')

par.paramvsRpart(Npx, pxlen, 15, 10, 30, 11,
                 mf.genPyramidTip, 65, 'angle', np.pi/8, np.pi/2, 1,
                 1, lambda surf: par.coverage(surf, thres), 'covvsR.dat')