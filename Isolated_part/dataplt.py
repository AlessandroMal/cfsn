import numpy as np
import matplotlib.pyplot as plt
from io import StringIO

def dataforplot(filename):
    inpt=open(filename)
    file=inpt.read()
    
    comment1=file.find('#')
    end1=file.find('\n', comment1+1)
    comment2=file.find('#', end1+1)
    end2=file.find('\n', comment2+1)
    if file.find('tip_r',comment2+1,end2)!=-1: tipparname=r'$R_{tip}/\langle R_{p}\rangle :$'
    else: tipparname=r'$\alpha :$' 
    return np.fromstring(file[:comment1-1], sep=' '), np.fromstring(file[end1+1:comment2-1], sep=' '), np.loadtxt(StringIO(file[end2+1:])), tipparname
#main=np.fromstring(file[:ind1-1], sep=' ')
#tip_par=np.fromstring(file[ind1+1:ind2-1], sep=' ')
#datas=np.loadtxt(StringIO(file[ind2+1:]))


def plotfile(f,rel=True):
    main,tip,data,name=dataforplot(f)
    if rel: 
        for i in range(len(tip)): data[i+2]=data[i+2]/data[1]

    if name.find('alpha')!=-1: lab=np.round(tip*180/np.pi, 1)
    else: lab=np.round(tip/main[2], 2)
    
    plt.figure()
    
    if f.find('vsN')!=-1:
        for i in range(len(tip)+1):
            if i==0 and not(rel): plt.plot(data[0]/main[4],data[1+i], label='surf')
            if i!=0: plt.plot(data[0]/main[4],data[1+i], label=name+str(lab[i-1]))
    
        plt.title(r'$Npx:$'+str(main[0])+',  ' + r'$\langle R_{part} \rangle /L_{px}=$'+str(round(main[2]/main[1] ,3))+',  '+r'$h_{tip}/L_{px}=$'+str(round(main[3]/main[1], 3)) )
        plt.xlabel(r'$N_{part}/N_{cp}$')
        plt.ylabel(f[f.find('/')+1 : f.find('vs')])
    
    else:
        for i in range(len(tip)+1):
            if i==0 and not(rel): plt.plot(data[0]/main[1],data[1+i], label='surf')
            else: plt.plot(data[0]/main[1],data[1+i], label=name+str(lab[i-1]))
    
        plt.title(r'$Npx:$'+str(main[0])+',  ' + r'$N_{part}=$'+str(main[2])+',  '+r'$h_{tip}/L_{px}=$'+str(round(main[3]/main[1], 3)) )
        plt.xlabel(r'$R_{part}/L_{px}$')
        plt.ylabel(f[f.find('/')+1 : f.find('vs')])
    
    
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()
    

plotfile('cone_tip/covvsR.dat', False)