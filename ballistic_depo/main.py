import mapsfunction as mf
#import parameters as par
#import scipy.ndimage.morphology as mph
import numpy as np
#import matplotlib.pyplot as plt
import part_num

#plt.close('all')  # chiude tutte le figure aperte

Npx = 800  # resolution of the map
pxlen = 1  # nm, phyisical length of a pixel in the map
#thres = 0  # nm
#Npart = 100  # number of particles on the map
#h = 80
#r = 4

R_mean_real = 30
R_std_real = 10
#R_mu_real = np.log(R_mean_real / np.sqrt(1 + R_std_real **
#                                         2 / R_mean_real**2))  # recalculated gaussian
#R_sigma_real = np.sqrt(np.log(1 + (R_std_real/R_mean_real)**2)) # reculaculated gaussian

# MAP-------------------------------------------------

#z = mf.genFlat(Npx)
#z = mf.genIsolLogNormSph(z, pxlen, Npart, R_mean_real, R_std_real)
#z = mf.genLogNormSph(z, pxlen, Npart, R_mean_real, R_std_real)
#z,R = mf.genLogNormSolidSph(z,pxlen,Npart,R_mean_real, R_std_real)

# Tip------------------------------------------------
# occhio che h/a>>pxlen

#tip = mf.genParabolicTip(pxlen, h, r=r)
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

part_num.NpartDep(Npx, pxlen, 3, 5, 500, 4, R_mean_real, R_std_real, 'N')
#part_num.VfracDep(Npx, pxlen, 15, 5, 10000, 30, R_mean_real, R_std_real)
#part_num.GrowthDep(Npx, pxlen, 15, 5, 10000, 30, R_mean_real, R_std_real)
# part_num.plotTipDep(z, pxlen, h, R_mu_real, R_sigma_real, Npart)

#mf.plotmap_fromfile('maps/lognorm_800_1_2.dat')
#part_num.plotpar_fromfile('N_relvsN_real.dat', r'$N_{real}$', r'$N_{est} / N_{real}$', xscale='log', yscale='log')