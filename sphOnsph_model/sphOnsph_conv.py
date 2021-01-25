import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
#import time

plt.close('all') #chiude tutte le figure aperte
#probabilmente devo avere Npxtip<280 circa (16gb ram con 25% circa occupato)
#278 funziona, 296 no
Npx=3000 #resolution of the map
pxlen=1 #nm, phyisical length of a pixel in the map
thres=0 #nm
s=1 #scaling factor to eventually avoid memory error (dilation)

rtip=36
r_part_min=s*12
r_part_max=s*108

htip=r_part_max*2*1.02

h_arr =np.array([])
r_arr =np.array([])
A_arr =np.array([])
V_arr =np.array([])
dV_arr=np.array([])

posmax=np.array([])
profiles=[]
height=[]

dist=2*(rtip + r_part_max)
dist*=np.sqrt(2)*1.0

z=mf.genFlat(Npx)
z=mf.genHexCap(z,pxlen,dist,r_part_min,r_part_max, elem='sphere',
               xmin=dist, xmax=len(z)*pxlen-dist,
               ymin=dist, ymax=len(z)*pxlen-dist)
        
tip=mf.genSemisphTip(pxlen,htip,r=rtip)
img = mph.grey_dilation(z, structure=-tip)
mf.plotfalsecol(z,pxlen)
mf.plotfalsecol(img,pxlen)
mf.plotfalsecol(tip,pxlen)
if 2*r_part_max>htip: print('Warning: possible spikes higher than semisphere tip')
obj = mf.identObj(img,thres)
for o in obj:
#       mf.plotfalsecol(o,pxlen)
    h,r,A,V,e=par.capPar(o,pxlen,thres)
    h_arr=np.append(h_arr, h)
    r_arr=np.append(r_arr, r)
    A_arr=np.append(A_arr, A)
    V_arr=np.append(V_arr, V)
    dV_arr=np.append(dV_arr, 1 - 6*V/(np.pi*h**3 + 3*h*A))
        

np.savetxt('sphOnsph.dat', np.array([h_arr, r_arr, A_arr, V_arr, dV_arr]),
           header='R_tip='+str(rtip)+'; Npx='+str(Npx)+'; pxlen='+str(pxlen)+'\n on rows: h, r, A, V, dV')
