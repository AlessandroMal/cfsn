import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
#import time

plt.close('all') #chiude tutte le figure aperte
#probabilmente devo avere Npxtip<280 circa (16gb ram con 25% circa occupato)
#278 funziona, 296 no
Npx=1500 #resolution of the map
pxlen=1 #nm, phyisical length of a pixel in the map
thres=0 #nm
s=1 #scaling factor to eventually avoid memory error (dilation)

rtip=200
r_part_min=s*10
r_part_max=s*20

htip=r_part_max*1.02

h_arr =np.array([])
r_arr =np.array([])
A_arr =np.array([])
V_arr =np.array([])
dV_arr=np.array([])

posmax=np.array([])
profiles=[]
height=[]

dist=2*np.sqrt(r_part_max**2 + 2*rtip*r_part_max)
#dist=2*(rtip + r_part_max)
#dist*=1.03

z=mf.genFlat(Npx)
z=mf.genHexCap(z,pxlen,dist,r_part_min,r_part_max,
               xmin=dist/2, xmax=len(z)*pxlen-dist/2,
               ymin=dist/2, ymax=len(z)*pxlen-dist/2)

mf.plotfalsecol(z,pxlen)

tip=mf.genSemisphTip(pxlen,htip,r=rtip)
img = mph.grey_dilation(z, structure=-tip)
mf.plotfalsecol(img,pxlen)
mf.plotfalsecol(tip,pxlen)
if r_part_max>htip: print('Warning: possible spikes higher than semisphere tip')
obj = mf.identObj(img,thres)[0]
for o in obj:
#    mf.plotfalsecol(o,pxlen)
    h,r,A,V,e=par.capPar(o,pxlen,thres)
    h_arr=np.append(h_arr, h)
    r_arr=np.append(r_arr, r)
    A_arr=np.append(A_arr, A)
    V_arr=np.append(V_arr, V)
    dV_arr=np.append(dV_arr, 1 - 6*V/(np.pi*h**3 + 3*h*A))
        

np.savetxt('capOnsph5.dat', np.array([h_arr, r_arr, A_arr, V_arr, dV_arr]),
           header='R_tip='+str(rtip)+'; Npx='+str(Npx)+'; pxlen='+str(pxlen)+'\n on rows: h, r, A, V, dV')
