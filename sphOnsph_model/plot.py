import numpy as np
import matplotlib.pyplot as plt

inp=open('dilation_measures.dat') #leggo i dati

lin=inp.readline()
R_rel=np.fromstring(lin[lin.find(':')+2 : ])
lin=inp.readline()
pxlen=np.fromstring(lin[lin.find(':')+2 : ])
lin=inp.readline()
V_dig=np.fromstring(lin[lin.find(':')+2 : ])
lin=inp.readline()
V_dil=np.fromstring(lin[lin.find(':')+2 : ])
lin=inp.readline()

R=R_rel[-1]
R_rel=np.delete(R_rel, -1) / R

x_prof= []
z_prof= []
lin=inp.readline()
x_prof.append(np.fromstring(lin))
lin=inp.readline()
z_prof.append(np.fromstring(lin))
lin=inp.readline()
for i in range(len(R_rel)):
    lin=inp.readline()
    x_prof.append(np.fromstring(lin))
    lin=inp.readline()
    z_prof.append(np.fromstring(lin))
    inp.readline()
inp.close()

plt.figure()
plt.title(r'Volume ratios, $d/L_{px}=$'+str(round(2*R/pxlen[-1], 3)))
plt.plot(R_rel,V_dig[len(pxlen)-1::len(pxlen)], label=r'$V_{map}$')
plt.plot(R_rel,V_dil[len(pxlen)-1::len(pxlen)], label=r'$V_{dil}$')
plt.xlabel(r'$ R_{tip} / R_{part} $')
plt.xlabel(r'$ V_{tip} / V_{real} $')
plt.legend()
plt.grid()
plt.show()


plt.figure()
plt.title(r'particle profile at $d/L_{px}=$'+str(round(2*R/pxlen[-1], 3)))
for i in range(0,len(z_prof)):
    if i==0: plt.plot(x_prof[i], z_prof[i],  label='digital map')
    else: plt.plot(x_prof[i], z_prof[i], label=r'$ R_{tip} / R_{part} = $'+str(round(R_rel[i], 2)))

plt.xlabel('r [npx]')
plt.xlabel('z [npx]')
plt.legend()
plt.grid()
plt.show()