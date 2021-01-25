import numpy as np
import matplotlib.pyplot as plt


h,r,A,V,dV=np.loadtxt('sphOnsph.dat')

inp=open('sphOnsph.dat') #leggo i dati
file=inp.readline()
inp.close()

start=file.find('R_tip=') +6
end  =file.find(';', start)
Rtip=float(file[start:end])
start=file.find('Npx=', end) +4
end  =file.find(';', start)
Npx=float(file[start:end])
start=file.find('pxlen=', end) +6
end  =file.find(';', start)
pxlen=float(file[start:end])

plt.figure()
plt.title(r'$V_{cap}$ error function')
plt.scatter(Rtip/h,dV)
plt.xlabel(r'$ R_{tip} / R_{part} $')
plt.ylabel(r'$ \Delta V $')
plt.grid()
plt.show()

plt.figure()
plt.title('Deformation')
plt.scatter(Rtip/h,r/h)
plt.xlabel(r'$ R_{tip} / R_{part} $')
plt.ylabel(r'$ R_{dil} / R_{true} $')
plt.grid()
plt.show()