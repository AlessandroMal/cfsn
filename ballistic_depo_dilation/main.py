import mapsfunction as mf
#import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
import part_num
import random

plt.close('all')  # chiude tutte le figure aperte

Npx = 1000  # resolution of the map
pxlen = 0.5  # nm, phyisical length of a pixel in the map
#thres = 0  # nm
#Npart = 100000  # number of particles on the map
Npart = 100000
#h = 80
#r = 4

h = 20
#R_tip = np.linspace(1, 25, 10)
R_tip = np.linspace(0.01, 3, 10)

N_sample = 10

R_median_real = 5
R_std_real = 1
R_mu_real = np.log(R_median_real)
R_sigma_real = np.sqrt(np.log(1/2 + np.sqrt(1/4 + np.exp(-2*R_mu_real)*R_std_real**2)))

# MAP-------------------------------------------------

#z = mf.genFlat(Npx)
#z = mf.genIsolLogNormSph(z, pxlen, Npart, R_mean_real, R_std_real)
#z = mf.genLogNormSph(z, pxlen, Npart, R_mean_real, R_std_real)
#z, R = mf.genLogNormSolidSph(z,pxlen,Npart,R_mu_real, R_sigma_real)

# Tip------------------------------------------------
# occhio che h/a>>pxlen

#tip = mf.genParabolicTip(pxlen, h, r=5)
# tip=mf.genPyramidTip(pxlen,50,angle=np.pi/4)
# tip=mf.genConeTip(pxlen,50,angle=np.pi/2)
# tip=mf.genSemisphTip(pxlen,50,r=20)

# IMG------------------------------------------------

#img = mph.grey_dilation(z, structure=-tip)

# PLOT-----------------------------------------------

#mf.plotview(z, pxlen, 30, 30)
# mf.plotview(z,pxlen,90,0)
# mf.plotview(z,pxlen,0,90)
#mf.plotfalsecol(z, pxlen)
#mf.plotview(img, pxlen, 30, 30)
#mf.plotview(img,pxlen,90,0)
#mf.plotview(img,pxlen,0,90)
#mf.plotfalsecol(img, pxlen)
#mf.plotview(tip, pxlen, 30, 30)
# mf.plotview(tip,pxlen,90,0)
# mf.plotview(tip,pxlen,0,90)

# PARTNUM--------------------------------------------

# R_mean, R_std, R_mu, R_sigma = part_num.reconstLogNorm(
#     z, pxlen, thres, Npart, R_mu_real, R_sigma_real)
# print('R_mean_real:', R_mean_real, 'R_std_real:', R_std_real,
#       '\nR_mu_real:', R_mu_real, 'R_sigma_real:', R_sigma_real)
# print('R_mean:', R_mean, 'R_std:', R_std, '\nR_mu:', R_mu, 'R_sigma:', R_sigma)
#print(part_num.partNum(z, pxlen, R_mu_real, R_sigma_real))

#part_num.partDep(Npx, pxlen, 3, 5, 500, 4, R_mean_real, R_std_real, 'N')
#part_num.VfracDep(Npx, pxlen, 15, 5, 10000, 30, R_mean_real, R_std_real)
#part_num.GrowthDep(Npx, pxlen, 15, 5, 10000, 30, R_mean_real, R_std_real)
#part_num.plotTipDep(z, pxlen, h, R_mu_real, R_sigma_real, Npart)


part_num.saveImg(Npx, pxlen, Npart, R_median_real, R_std_real, R_mu_real, R_sigma_real, R_tip, h, N_sample, note='zoom')
part_num.RMSTipDep(Npart, R_median_real, R_std_real, R_tip, N_sample, note='zoom')
part_num.specAreaTipDep(pxlen, Npart, R_median_real, R_std_real, R_tip, N_sample, note='zoom')
plt.show()

#mf.plotmap_fromfile('maps/lognorm_800_1_2.dat')
#part_num.plotpar_fromfile('N_relvsN_real.dat', r'$N_{real}$', r'$N_{est} / N_{real}$', xscale='log', yscale='log')
