import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
from scipy.ndimage import median_filter
import numpy as np
import matplotlib.pyplot as plt
import part_num
import os
import dill
import random

np.random.seed(79223)
random.seed(79223)

plt.close('all')  # chiude tutte le figure aperte

Npx = 1000  # resolution of the map
pxlen = 0.25  # nm, phyisical length of a pixel in the map
thres = 0 # nm
Npart = 100 # number of particles on the map
h = 15
r = 5

sigma_noise = 0

R_median_real = 5
R_std_real = 1
R_mu_real = np.log(R_median_real)
R_sigma_real = np.sqrt(np.log(1/2 + np.sqrt(1/4 + np.exp(-2*R_mu_real)*R_std_real**2)))

# MAP-------------------------------------------------

z = mf.genFlat(Npx)
z = mf.genIsolLogNormSph(z, pxlen, Npart, R_mu_real, R_sigma_real)
z = mf.genNoise(z, Npx, 0, sigma_noise)
#z = mf.genLogNormSph(z, pxlen, Npart, R_mean_real, R_std_real)
#z,R = mf.genLogNormSolidSph(z,pxlen,Npart,R_mean_real, R_std_real)

#z = median_filter(z, 3)

# Tip------------------------------------------------
# occhio che h/a>>pxlen

tip = mf.genParabolicTip(pxlen, h, r=r)
# tip=mf.genPyramidTip(pxlen,50,angle=np.pi/4)
# tip=mf.genConeTip(pxlen,50,angle=np.pi/2)
# tip=mf.genSemisphTip(pxlen,50,r=20)

# IMG------------------------------------------------

img = mph.grey_dilation(z, structure=-tip)

# PLOT-----------------------------------------------

mf.plotview(z, pxlen, 30, 30)
# mf.plotview(z,pxlen,90,0)
# mf.plotview(z,pxlen,0,90)
mf.plotfalsecol(z, pxlen)
mf.plotview(img, pxlen, 30, 30)
#mf.plotview(img,pxlen,90,0)
#mf.plotview(img,pxlen,0,90)
mf.plotfalsecol(img, pxlen)
#mf.plotview(tip, pxlen, 30, 30)
# mf.plotview(tip,pxlen,90,0)
# mf.plotview(tip,pxlen,0,90)

# PARTNUM--------------------------------------------

R_median, R_mu, R_sigma = part_num.reconstLogNorm(
    z, pxlen, r, sigma_noise, thres, Npart, R_median_real, R_mu_real, R_sigma_real)

# print('R_mean_real:', R_mean_real, 'R_std_real:', R_std_real,
#       '\nR_mu_real:', R_mu_real, 'R_sigma_real:', R_sigma_real)
# print('R_mean:', R_mean, 'R_std:', R_std, '\nR_mu:', R_mu, 'R_sigma:', R_sigma)
# print(part_num.partNum(z, pxlen, R_mu_real, R_sigma_real))

# part_num.plotNpartDep(Npx, pxlen, 30000, 10000, R_mean_real, R_std_real, R_mu_real, R_sigma_real)
# part_num.plotVfracDep(Npx, pxlen, 50000, 1000, R_mean_real, R_std_real, R_mu_real, R_sigma_real)
# part_num.plotGrowthDep(Npx, pxlen, 50000, 1000, R_mean_real, R_std_real, R_mu_real, R_sigma_real)
# part_num.plotTipDep(z, pxlen, h, R_mu_real, R_sigma_real, Npart)

# SAVE-----------------------------------------------

os.makedirs('results', exist_ok=True)
dill.dump_session('results/saved_session')
#dill.load_session('results/saved_session')

figs = [plt.figure(n) for n in plt.get_fignums()]
for i, fig in enumerate(figs):
    fig.savefig('results/fig'+str(i), format='png')