import mapsfunction as mf
import parameters as par
#import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
#from scipy.ndimage import gaussian_filter
#import time

plt.close('all') #chiude tutte le figure aperte

def calcAnsph(r,L):
    av=5/3*np.pi*r**3/L**2
    stdev=np.sqrt(np.pi/3*r**4/L**2*(17/2-25/3*np.pi*r**2/L**2))
    sk=av/stdev**3*(213/70*r**2 - 3*stdev**2 - av**2)
    params = {'mean': av,
              'std': stdev,
              'skew': sk,
              'V': 5/3*np.pi* r**3,
              'specArea': 1+ 3*np.pi* r**2/L**2,
              'coverage': np.pi* (r/L)**2}
    return params


res=[]
numeric=[]
analytic=[]
thres=0
L=100
q=4/10
r=q*L
for i in range(10,501,20): #iterazione su risoluzione
    pxlen=L/i
    res.append(i)
    
    z=mf.genFlat(i)
    z=mf.genSphere(z,pxlen,np.array([L/2,L/2]),np.array([r]))
    numeric.append(par.calcParams(z,pxlen,thres))
    analytic.append(calcAnsph(r,L))
        
for i in numeric[0]:
    ntemp = []
    atemp = []
    for j in range(len(numeric)):
        ntemp.append(numeric[j][i])
        atemp.append(analytic[j][i])
        
    plt.figure()
    plt.title('Analytic-measured paramteres comparison for single sphercal particle')
    plt.plot(res, ntemp, marker='.', color='green', label='measured')
    plt.plot(res, atemp, color='red', label='analytic')   
    plt.xlabel('Resolution of L (d='+str(q*2)+'L)')
    plt.ylabel(i)
    plt.legend()
    plt.grid()
    print(i)
    print('at lower res (analytic-measured) : ',analytic[0][i],' ',numeric[0][i],' ',(numeric[0][i]-analytic[0][i])/numeric[0][i])
    print('at higher res (analytic-measured): ',analytic[-1][i],' ',numeric[-1][i],' ',(numeric[-1][i]-analytic[-1][i])/numeric[-1][i])

#TIME---------------------------------------------
'''
v=[]
a=100
b=1000
step=100
b+=1
for i in range(a,b,step):
    z=mf.genFlat(i)
    start=time.perf_counter()
    par.specArea(z,pxlen)
    v.append(time.perf_counter()-start)
    
plt.figure(figsize=(20,20))
plt.plot(np.arange(a,b,step), v, marker='.')
plt.plot(np.arange(a,b,step),6.5/9/10000*np.arange(a,b,step)**2)
plt.grid()
'''  