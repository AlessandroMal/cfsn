import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
#from scipy.ndimage import gaussian_filter
#import time

plt.close('all') #chiude tutte le figure aperte
#probabilmente devo avere Npxtip<280 circa (16gb ram con 25% circa occupato)
#278 funziona, 296 no

def volsphere(R): return 5/3*np.pi*R**3
def volsemisph(R): return 2/3*np.pi*R**3

L=100 #lato mappa
q=1/10 #rapporto raggio/lato
Npx=np.linspace(100,500,8)
#Npx=np.linspace(100,42,2)

pxlen=L/Npx
thres=0 #soglia
#-----------------------------------------------
R=q*L
htip=2*R*1.02
R_tip=np.linspace(R/5,3*R, 8)
#R_tip=np.linspace(3*R,3.1*R, 2)

V_red_calc= []
V_red_dil = []
height = []
profiles = []
posmax = []
for rtip in R_tip:
    for i in range(len(Npx)): #iterazione su risoluzione
        print('R_tip=', rtip , ' Npx=', int(Npx[i]))
        z=mf.genFlat(int(Npx[i]))
        z=mf.genSphere(z,pxlen[i],np.array([Npx[i]*pxlen[i]/2, Npx[i]*pxlen[i]/2]),np.array([R]))
        obj = mf.identObj(z,thres)[0]
        V_red_calc.append(par.V(obj, pxlen) / volsphere(R))
        
        if i==len(Npx)-1 and rtip==R_tip[0]:
            maxh=0
            posmax.append(0)
            for x in range(np.shape(obj)[0]):
                height.append(np.amax(obj[0:,x])/R)
                if maxh<height[-1]:
                    maxh=height[-1]
                    posmax[-1]=x
            profiles.append(np.array(height))
            height.clear()
        
        tip=mf.genSemisphTip(pxlen[i],htip,r=rtip)
        z = mph.grey_dilation(z, structure=-tip)
        mf.plotfalsecol(z,pxlen[i])
        obj = mf.identObj(z,thres)[0]
        V_red_dil.append(par.V(obj, pxlen[i]) / volsphere(R))
        
        if i==len(Npx)-1:
            maxh=0
            posmax.append(0)
            for x in range(np.shape(obj)[0]):
                height.append(np.amax(obj[0:,x])/R)
                if maxh<height[-1]:
                    maxh=height[-1]
                    posmax[-1]=x
            profiles.append(np.array(height))
            height.clear()

R_tip=np.append(R_tip, R)
print('printing data in dilation_measures.dat ...')
out=open('dilation_measures.dat', 'w')
out.write('r_tip (R_part in last pos): '+str(R_tip)+'\n')
out.write('Lpx: '+str(pxlen)+'\n')
out.write('digital map volume: '+str(np.array(V_red_calc))+'\n')
out.write('dilation map volume: '+str(np.array(V_red_dil))+'\n')
out.write('x_surf, surface profile, x dilated, map profile\n')
out.write(str(pxlen[-1] * (np.arange(0,len(profiles[0])) - posmax[0]) )+'\n')
out.write(str(profiles[0])+'\n'+'\n')
for i in range(1,len(R_tip)+1):
    out.write(str(pxlen[-1] * (np.arange(0,len(profiles[i])) - posmax[i]) )+'\n')
    out.write(str(profiles[i])+'\n'+'\n')
out.close()
