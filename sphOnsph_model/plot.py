import numpy as np
import matplotlib.pyplot as plt

inp=open('conv_sphere.dat') #leggo i dati
file=inp.read()
inp.close()

start=file.find('[') +1
end  =file.find(']', start)
R_rel=np.fromstring(file[start:end].replace('\n',' '), sep=' ')
start=file.find('[', end+1) +1
end  =file.find(']', start)
pxlen=np.fromstring(file[start:end].replace('\n',' '), sep=' ')
start=file.find('[', end+1) +1
end  =file.find(']', start)
V_dig=np.fromstring(file[start:end].replace('\n',' '), sep=' ')
start=file.find('[', end+1) +1
end  =file.find(']', start)
V_dil=np.fromstring(file[start:end].replace('\n',' '), sep=' ')

R=R_rel[-1]
R_rel=np.delete(R_rel, -1) / R

x_prof= []
z_prof= []
for i in range(len(R_rel)+1):
    start=file.find('[', end+1) +1
    end  =file.find(']', start)
    x_prof.append(np.fromstring(file[start:end].replace('\n',' '), sep=' '))
    start=file.find('[', end+1) +1
    end  =file.find(']', start)
    z_prof.append(np.fromstring(file[start:end].replace('\n',' '), sep=' '))

plt.figure()
plt.title(r'Volume ratios, $d/L_{px}=$'+str(round(2*R/pxlen[-1], 3)))
#plt.plot(R_rel,V_dig* R_rel**0, label=r'$V_{map}$')
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
    else: plt.plot(x_prof[i], z_prof[i], label=r'$ R_{tip} / R_{part} = $'+str(round(R_rel[i-1], 2)))

plt.xlabel('r [npx]')
plt.xlabel('z [npx]')
plt.legend()
plt.grid()
plt.show()