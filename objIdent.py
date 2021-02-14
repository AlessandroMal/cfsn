def identObj(z, thres):
    z_binary = (z>thres) #vero se elemento e piu grande della soglia
    z_labeled = scipy.ndimage.label(z_binary)[0] #numera le singole particelle
    obj_ind = scipy.ndimage.find_objects(z_labeled) #trova gli indici delle particelle
    
    obj_list = []
    
    # prevent overlap of single structures
    for i in range(len(obj_ind)):
        z_single_obj = z.copy()
        z_single_obj[np.where(z_labeled!=i+1)] = 0
        
        # remove zero rows and columns
        z_single_obj = z_single_obj[~np.all(z_single_obj == 0, axis=1)]
        z_single_obj = z_single_obj[:, ~np.all(z_single_obj == 0, axis=0)]
        
        obj_list.append(z_single_obj)

    return obj_list, z_labeled, obj_ind
