import numpy as np
import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph

def reconstLogNorm(z, pxlen, thres, N_part):
    z_obj_list, z_labeled = mf.identObj(z, thres)
    mf.plotThres(z, z_labeled, pxlen, 'found ' + str(len(z_obj_list)) + ' of ' + str(N_part) + ' particles, thres=' + str(thres)+'nm')
    
    R_list = [(np.max(np.shape(obj_i))-1)*pxlen/2 for obj_i in z_obj_list]
    
    R_mean = np.mean(R_list)
    R_std = np.std(R_list, ddof=1)
    
    R_mu = np.log(R_mean / np.sqrt(1 + R_std**2 / R_mean**2)) # recalculated gaussian
    R_sigma = np.sqrt(np.log(1 + (R_std/R_mean)**2)) # reculaculated gaussian

    return R_mean, R_std, R_mu, R_sigma
    
def partNum(z, pxlen, R_mu, R_sigma):
    A = np.len(z)**2 * pxlen**2
    V = par.V(z, pxlen)
    
    V_mean = 4/3 * np.pi * np.exp(3*(2*R_mu + k*R_sigma**2)/2)