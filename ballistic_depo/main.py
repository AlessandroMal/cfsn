import mapsfunction as mf
import parameters as par
#import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
import part_num

#plt.close('all')  # chiude tutte le figure aperte

Npx = 128  # resolution of the map
pxlen = 1  # nm, phyisical length of a pixel in the map
#thres = 0  # nm
#Npart = 1  # number of particles on the map
#h = 80
#r = 4

R_median = 7
R_std = 0.33
R_mu_real = np.log(R_median)

#R_sigma_real = np.sqrt(np.log(1/2 + np.sqrt(1/4 + np.exp(-2*R_mu_real)*R_std**2)))
#beamw= 2*10**7
#print(R_sigma_real)
#R_mu_real = np.log(R_mean_real / np.sqrt(1 + R_std_real **
#                                         2 / R_mean_real**2))  # recalculated gaussian
#R_sigma_real = np.sqrt(np.log(1 + (R_std_real/R_mean_real)**2)) # reculaculated gaussian

# MAP-------------------------------------------------

#z = mf.genFlat(Npx)
#z = mf.genIsolLogNormSph(z, pxlen, Npart, R_mean_real, R_std_real)
#z = mf.genLogNormSph(z, pxlen, Npart, R_mean_real, R_std_real)
#z= mf.genUnifSolidSph(z,pxlen,Npart,30,30, xmin=0, xmax=140, ymin=60, ymax=110)
#z= mf.genUnifSolidSph(z,pxlen,Npart,30,30, xmin=0, xmax=140, ymin=60, ymax=110)
#z= mf.genUnifSolidSph(z,pxlen,Npart,30,30, xmin=50, xmax=140, ymin=60, ymax=110)
#print('top', np.amax(z))
#np.savetxt('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(Npart)[:len(str(Npart))-2]+'.dat', z, header='V='+str(np.sum(4/3 * np.pi * R_l**3))+'; Npx,pxlen,Npart in filename')
#xcorr, Ccorr= par.G_profile(z,2,250)
#xcor2,Gcorr=par.G_profile(z,1,600)
#plt.figure()
#plt.xscale('linear')
#plt.yscale('log')
#plt.plot(xcorr,Ccorr)
#plt.grid()
#plt.show()

#plt.figure()
#plt.xscale('log')
#plt.yscale('log')
#plt.plot(xcor2,Gcorr)
#plt.grid()
#plt.show()

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
#mf.plotview(z,pxlen,90,0)
#mf.plotview(z,pxlen,0,90)
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

#out=open('parametri.mat', 'w')
#out.write('2^16 1 1024 1024 1024 1024')
#out.close()

#part_num.partDep(Npx, pxlen, 1, 300, 30000, 4, R_mu_real, R_std, savefile=True)
part_num.partDep(Npx, pxlen, 10, 1000, 20000000, 60, R_mu_real, R_std)
#part_num.plotTipDep(z, pxlen, h, R_mu_real, R_sigma_real, Npart)

#mf.plotmap_fromfile('maps/lognorm_1024_1_300.dat')
#mf.plotmap_fromfile('maps/lognorm_1024_1_3000.dat')
#mf.plotmap_fromfile('maps/lognorm_1024_1_30000.dat')
#mf.plotmap_fromfile('maps/lognorm_1024_1_300000.dat')
#part_num.plotpar_fromfile('N_relvsN_real.dat', r'$N_{real}$', r'$N_{est} / N_{real}$', xscale='log', yscale='linear', a=80,b=5000)
#part_num.plotpar_fromfile('V_relvsN_real.dat', r'$N_{real}$', r'$V_{est} / V_{real}$', xscale='log', yscale='linear', a=80,b=5000)
#part_num.plotpar_fromfile('rmsvsN_real.dat', r'$N_{real}$', r'$RMS [L_{px}]$', xscale='log', yscale='log', a=10,b=1000)
