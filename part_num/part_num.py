import numpy as np
import mapsfunction as mf
import scipy.ndimage.morphology as mph

def reconstLogNorm(z, pxlen, thres, N_part):
    z_obj_list, z_labeled = mf.identObj(z, thres)
    mf.plotThres(z, z_labeled, pxlen, 'found ' + str(len(z_obj_list)) + ' of ' + str(N_part) + ' particles, thres=' + str(thres)+'nm')
    
    R_list = [(np.max(np.shape(obj_i))-1)*pxlen/2 for obj_i in z_obj_list]
    
    R_mean = np.mean(R_list)
    R_std = np.std(R_list, ddof=1)

    print(R_mean, R_std)