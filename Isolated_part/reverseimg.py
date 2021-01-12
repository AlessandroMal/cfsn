import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from scipy.ndimage import gaussian_filter
#import time

plt.close('all') #chiude tutte le figure aperte
#probabilmente devo avere Npxtip<280 circa (16gb ram con 25% circa occupato)
#278 funziona, 296 no
Npx=2030 #resolution of the map
pxlen=5 #nm, phyisical length of a pixel in the map
thres=1 #nm
s=1 #scaling factor to eventually avoid memory error (dilation)

rspikemin=s*55
rspikemax=s*80
hspikemin=s*110
hspikemax=s*160

rtip=s*1500

htip=hspikemax*1.02
spikedist=2*np.sqrt(rtip**2 - (rtip-hspikemax)**2)
#MAP-------------------------------------------------
'''
z=mf.genFlat(Npx)
z=mf.genHexSpikes(z,pxlen,hspikemin,hspikemax,spikedist,rspikemin,rspikemax,0.8,0,
                  par='r', xmin=spikedist, xmax=len(z)*pxlen-spikedist,
                  ymin=spikedist, ymax=len(z)*pxlen-spikedist)

#Tip------------------------------------------------
#occhio che h/a>>pxlen
tip=mf.genSemisphTip(pxlen,htip,r=rtip)

#IMG------------------------------------------------

img = mph.grey_dilation(z, structure=-tip)

#PLOT-----------------------------------------------

mf.plotfalsecol(z,pxlen)
#mf.plotview(z,pxlen,90,-90)

mf.plotfalsecol(img,pxlen)
mf.plotfalsecol(tip,pxlen)

#---------------------------------------------------
if hspikemax>rtip: print('Warning: possible spikes higher than semisphere tip')


obj = mf.identObj(img,thres)
#for i in obj: mf.plotfalsecol(i,pxlen)
filtobj= par.revimg_filter(obj,pxlen,thres,0.3,1.5,0.2)
#for i in obj: print(par.capPar(i,pxlen,thres))
#for i in filtobj: mf.plotfalsecol(i,pxlen)

h=[]
V=[]
A=[]
for i in filtobj:
    h.append(par.capPar(i,pxlen,thres)[0]) #h è 0, A è 2, V è 3
    V.append(par.capPar(i,pxlen,thres)[3])
    A.append(par.capPar(i,pxlen,thres)[2])

np.savetxt('revimg.dat', (np.array(h), np.array(V), np.array(A)),
           header='h,V,A; Npx='+str(Npx)+'; Rtip='+str(round(rtip,1)))
'''
def V_h(h_obj, R): return np.pi/3*h_obj**2 *(3*R-h_obj)
def A_h(h_obj, R): return np.pi*h_obj *(2*R-h_obj)
def V_A(A_obj, R): return 2/3*np.pi*R**3 *(1 - np.sqrt(1-A_obj/R**2 /np.pi) - A_obj/R**2 /np.pi/2) *(2 + np.sqrt(1-A_obj/R**2 /np.pi))

def regression(f, filename):
    parh,parv,para= np.loadtxt(filename)
    if str(f)[10]=='V': y=parv
    if str(f)[10]=='A': y=para
    if str(f)[12]=='h': x=parh
    if str(f)[12]=='A': x=para
    
    if str(f)[10:13]=='V_A': R_opt, R_var= curve_fit(f, x, y, bounds=(np.sqrt(min(x)/np.pi),1.5*rtip) )
    else: R_opt, R_var= curve_fit(f, x, y)
    print(str(f)[10]+'vs '+str(f)[12]+' curve:')
    print('Rtip=',*R_opt)
    print('Rstd=',np.sqrt(*R_var[0]))
    plt.figure()
    plt.title(r'Regression optimizing $R_{tip}$')
    plt.scatter(x,y, label='measures from caps')
    xrgr=np.linspace(min(x),max(x),100)
    plt.plot(xrgr,f(xrgr,*R_opt),color='red', label='regression')
    plt.xlabel(str(f)[12])
    plt.ylabel(str(f)[10])
    plt.grid()
    plt.legend()
    plt.show()

regression(V_h, 'reverse_img/revimg.dat')
regression(A_h, 'reverse_img/revimg.dat')
regression(V_A, 'reverse_img/revimg.dat')
