# standard modules
import os
import pickle
import time

import numpy as np
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as scidist
from md_fcns import fill_pocket, plot_heatmap, write_mdpocket

from AS_fcn import combine_prot_lig, genCommunityPocket, runAlphaSpace
# out module
from DTM_old.md_classes import MDpocket
from retired.AlphaSpace import readPDB


def genCMDpocket(pdblist, 
                param_file = "AS_param.txt", option_file = "AS_option.txt",
                metric = 'jaccard', cmethod = 'distance',
                sieve = False,  sieveopt = "random", step = 20, 
                isScan = False, scanlist = [],
                isdump = True, isplot = False, isdebug = False, istime = False):
    """Function to run pocket for a list of pdbs from MD simulation with 
    Sieved option
    
    :param pdblist: a list of pdb name with location
    :param minsphere: minimum of spheres forming a pocket
    :param metric: distance calculation method, "euclidean" or "jaccard"
    :param cmethod: cut tree method, "maxclust" or "distance"
    :param sieve: Ture or False for sieve
    :param sieveopt: option for sieve method, "random" or "even"
    
    :param isplot: True or False of plot heatmap
    :param isdebug: True or False for debugging
    :param issubpocket: True or False for subpocket
    :returns: md_pockets folder and mdpocketlist
    
    .. note:: 
        **cmethod**: (1)maxclust - cut based on nubmer of cluster 
        (average of pockets); (2) distance - cut based on the height 
        of the tree (0.7)
        
    .. note::
        ftdict[(SnapshotID, PocketID)] = [AtomArray(np.array), 
        AtomCenter(list), PocketCenter(list), Volume(f), Score(f)]
    
    """
        
    if istime: print("\nTime will be printed:")
    
    # dict for snapshot and pocket
    ftdict = {} # ftdict for clustering, key is (ssid, cpid)
    aftdict = {} # ftdict for all pocket, key is (ssid, pid)
    surfdict = {} 
    snapidx = [] # keys in order of ftdict
    asnapidx = [] # keys in order of aftdict
    
    npocket_list = [] # number of pocket for each snapshot
    
    nsnap = len(pdblist)
    
    if isScan:
        numscan = len(scanlist)
        pdblist = pdblist + scanlist
    else:
        numscan = 0
        
    # total num of snapshots
    numsnap = len(pdblist)
    
    ## get first pdb (template) info for future use
    if len(pdblist[0]) == 2:
        tmppdblist = combine_prot_lig(pdblist[0])
    elif len(pdblist[0]) == 1:
        tmppdblist = [lines for lines in open(pdblist[0][0])]
    
    tmppdb = readPDB(tmppdblist,isTER = True)
    natom = len(tmppdb.prot)
    
    if istime: t = time.time()
    
    olddir = os.getcwd()
    os.system('mkdir AlphaSpace')
    os.chdir('AlphaSpace')

    lig_resid = []
    lig_resnm = []
    
    # commuinity information
    commlist = {i: [] for i in range(numsnap)}


    ## runAlphaSpace for each snapshot
    for i in range(numsnap):
        pdb_file = pdblist[i]
        
        if len(pdb_file) == 2:
            pdblinelist = combine_prot_lig(pdb_file)
        elif len(pdb_file) == 1:
            pdblinelist = [lines for lines in open(pdb_file[0])]
        
        os.system('mkdir pdb' + str(i))
        os.chdir('pdb' + str(i))
        pdb, ss = runAlphaSpace(pdblinelist, 
                                param_file = param_file, 
                                option_file = option_file, 
                                default = False, isWrite=True)
        os.chdir('..')
        os.system("tar -cvf pdb" + str(i) + ".tar --remove-files pdb" + str(i))
        
        ssp = ss.sub_pockets
        npocket_list.append(len(ssp))    
        x = fill_pocket(natom, ssp)
        surfdict[(i)] = pdb.asa_prot
        
        # vol_simp
        # p.vol_simp(ss.simplex_idx_list, ss.prot_pdb)
        # vol_exc
        # p.vol_exc(ss.simplex_idx_list,ss.prot_pdb,ss.asa_radius2)

        pktComm = genCommunityPocket(pdb, ss)
        commlist[i] = pktComm

        iscomm = True
        if iscomm:
            pktComm = genCommunityPocket(pdb, ss)
            commlist[i] = pktComm
            commout = commFilter(ssp, pktComm, x)
            for j,p in enumerate(commout):
                snapidx.append((i,j))
                ftdict[(i,j)] = p

        for j,p in enumerate(ssp):
            asnapidx.append((i,j))
            p.set_vol()
            p.set_vresid()
            aftdict[(i,j)] = [x[j], list(p.atm_cntr()), \
                    list(p.vert_cntr()), p.space, p.score, p.dict_resid_vol]
        
        if lig_resid == []:
            for i in ss.resid:
                if i not in lig_resid:
                    lig_resid.append(i)
        
        if istime: 
            print(i, time.time() - t)
            t = time.time()
    
    os.chdir(olddir)
    
    #print commlist
    #print snapidx
    #print asnapidx
    #for i in ftdict: print i 

    ## get sieved snapidx (i,j) where i is in sieve_index. 
    ## If false sieve_snapidx = snapidx
    if sieve:
        if sieveopt == "random":
            sieve_index = genRandom(numsnap, step)
        elif sieveopt == "even":
            sieve_index = genEven(numsnap, step)
        else:
            print("Wrong sieve method")
        sieve_snapidx = []
        for (i,j) in snapidx:
            if i in sieve_index:
                sieve_snapidx.append((i,j))
    elif isScan:
        sieve_snapidx = []
        tmpidx = list(range(nsnap))
        for (i,j) in snapidx:
            if i in tmpidx:
                sieve_snapidx.append((i,j))
    else:
        sieve_snapidx = snapidx

    ## get atom array for each pocket each snapshot
    for i,j in enumerate(sieve_snapidx):
        if i == 0:
            pocket_atom_array = [ftdict[j][0]]
        else:
            pocket_atom_array = np.append(pocket_atom_array, [ftdict[j][0]],0)
    
    # np.savetxt("AtomArray_dataforplot.dat", pocket_atom_array)
    
    ## Clustering       
    ## distance metric
    if metric == 'euclidean':
        print('Metric is jaccard')
        data_dist = scidist.pdist(pocket_atom_array)
    elif metric == 'jaccard':
        print('Metric is jaccard')
        data_dist = scidist.pdist(pocket_atom_array, "jaccard")
    else:
        print("Error: wrong clustering method")
    # print data_dist.shape

    ## clustering algorithm
    zmat = hier.linkage(data_dist, method="average")
    
    ## cut the tree
    if cmethod == 'distance':
        #clust_dist =  0.6*max(zmat[:,2])
        #for i in np.arange(0.65, 0.95, 0.02):
        #    clust_dist = i*max(zmat[:,2])
        #    cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
        #    print i, max(cluster)
        clust_dist =  0.75 * max(zmat[:,2])
        cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
    elif cmethod == 'maxclust':
        clust_num = int(np.average(npocket_list))   
        cluster = hier.fcluster(zmat,clust_num,criterion='maxclust')
    else:
        print("Error: wrong cutoff method")
        
    if isplot:
        # np.savetxt("AtomArray_dataforplot.dat",pocket_atom_array)
        plot_heatmap(pocket_atom_array, figname ='md_clust_heatmap.png',
                     metric=metric)
   
    cluster = list(cluster)
    nmdpocket = max(cluster)   # number of pocket by MD clustering
    mdpocketidx = [[] for i in range(nmdpocket)]  # empty list of pkt idx
    
    #print cluster
    #print mdpocketidx


    ## asign snapshot index to each cluster
    ## cluster ID strarts from 1 in `cluster'
    for i,j in zip(sieve_snapidx, cluster):
        mdpocketidx[j-1].append(i)
    
    ## if sieve is true, assign the remain pocket to existing pocket
    ## 1. Get the average array for each existing cluster 
    ## 2. Assign the remain pocket to existing cluster
    
    if sieve:
        avgArray = genAvgArray(pocket_atom_array, cluster, nmdpocket)
        ## cluster ID starts from 1
        for i in snapidx:
            if i not in sieve_snapidx:
                tmp = avgArray - ftdict[i][0]
                norm = np.sum(tmp**2, axis = 1)
                #tmpidx = norm.index(min(norm))
                print(norm.shape)
                tmpidx = np.argmin(norm, axis = 0)
                mdpocketidx[tmpidx].append(i)
    
    if isScan:
        ## if sieve is true, the dtm pkt is also assigned based on the clustering 
        ## snapshot not on the overall snapshot

        avgArray = genAvgArray(pocket_atom_array, cluster, nmdpocket)
        ## cluster ID starts from 1
        for i in snapidx:
            if i not in sieve_snapidx:
                tmp = avgArray - ftdict[i][0]
                norm = np.sum(tmp**2, axis = 1)
                tmpidx = np.argmin(norm, axis = 0)
                mdpocketidx[tmpidx].append(i)
            
    if isdebug: print(mdpocketidx)

    cmdpocketidx = [[] for i in range(nmdpocket)]

    for idx,i in enumerate(mdpocketidx):
        print(i)

        for j in i:
            sid, cpid = j
            _, a =  commlist[sid][cpid]
            print(a)
            # this is for final calculation
            # given the core pocket, add aux pocket to that cluster
            #
            for k in a[0] + a[1]:
                print(a[0], a[1])
                if (sid, k) not in cmdpocketidx[idx]:
                    cmdpocketidx[idx].append((sid,k))
    for i in cmdpocketidx: print(i)
    # !!! assign aftdict to ftdict
    ftdict = aftdict
    mdpocketidx = cmdpocketidx

    
    mdpocketlist = [] # empty list of dtm pkt instance

    for i in mdpocketidx:
        tmp = dict((k, ftdict[k]) for k in i)
        if isdebug: print(tmp, i, nsnap)
        mdpocketlist.append(MDpocket(tmp, i, nsnap, surfdict, numscan=numscan))

    # if isdebug: print mdpocketlist
    
    ## order the mdpocketlist by the order of the occprobability
    occplist = []
    for p in mdpocketlist:
        occplist.append(p.occprob)
    
    mdpocketlist = [mdpocketlist[i] for i in list(np.argsort(occplist)[::-1])]
    mdpocketidx = [mdpocketidx[i] for i in list(np.argsort(occplist)[::-1])]
    write_mdpocket(mdpocketlist, tmppdb, nsnap)

    ## wirte all files for coding 
    if isdump:
        pickle.dump(ftdict, open('ftdict.pkl', 'w'))
        pickle.dump(mdpocketidx, open('mdpocketidx.pkl', 'w'))
        pickle.dump(surfdict, open('surfdict.pkl', 'w'))
        pickle.dump(tmppdb, open('tmppdb.pkl', 'w'))      
        pickle.dump(mdpocketlist, open('mdpocketlist.pkl', 'w'))
        pickle.dump(lig_resid, open('lig_resid.pkl', 'w'))
    
    ## surface state clustering
    #sscluster = getSSC(mdpocketlist)
    
    ## PPI modulation 
    #ppimod(mdpocketlist, lig_resid)
    
    return mdpocketlist



def commFilter(ssp, pktComm, x):
    """
    Input: Taking list of pocket, pocket community list and pocket atom array list
    Output: pocket list
    """

    ncp =  len(pktComm)

    out = [[] for i in range(ncp)]
    for i,p in enumerate(pktComm):
        _, a = p
        #print a[0], a[1], a[2]
        for j in a[0]:
            if out[i] == []:
                #out[i] = [np.int64(x[j]), [ssp[j].atm_cntr()], [ssp[j].vert_cntr()], ssp[j].space, ssp[j].score, ssp[j].dict_resid_vol]
                out[i] = [np.int64(x[j])]
            else:
                out[i][0] = np.int64(out[i][0]) | np.int64(x[j])
                #out[i][1].append(ssp[j].atm_cntr())
                #out[i][2].append(ssp[j].vert_cntr())
                #out[i][3] += ssp[j].space
                #out[i][4] += ssp[j].score
                #out[i][5] = {}
        #print out[i][1]
        #out[i][1] = list(np.average(out[i][1], axis = 0))
        #out[i][2] = list(np.average(out[i][2], axis = 0))
        #print out[i][1]
 
    #print out
    return out


def processCommSnap(commitem, ftdict, surfdict, cidx):
    """
    return size of comm pocket
    """
    if commitem == []:
        return [0,0,0]
    
    
    out = [[] for i in range(len(commitem))]
    for idx, item in enumerate(commitem):
        _, a = item
        tmp = []
        for i in a[0] + a[1] :
            if tmp == []:
                tmp = [ftdict[(cidx,i)][0], ftdict[(cidx,i)][3], ftdict[(cidx,i)][4]]
            else:
                tmp[0] = np.int64(tmp[0]) | np.int64(ftdict[(cidx,i)][0])
                tmp[1] += ftdict[(cidx,i)][3]
                tmp[2] += ftdict[(cidx,i)][4]
        sa = np.sum(tmp[0] * surfdict[cidx])
        tmp.append(sa)
        out[idx] = tmp[1:]
    return out


def processCommSnapAggregate(commitem, ftdict, surfdict, cidx):
    """
    return size of comm pocket
    """
    if commitem == []:
        return [0,0,0]
    tmp = []
    plist = []

    for idx, item in enumerate(commitem):
        _, a = item
        for i in a[0] + a[1]:
            if i not in plist:
                plist.append(i)
    for i in plist:
        if tmp == []:
            tmp = [ftdict[(cidx,i)][0], ftdict[(cidx,i)][3], ftdict[(cidx,i)][4]]
        else:
            tmp[0] = np.int64(tmp[0]) | np.int64(ftdict[(cidx,i)][0])
            tmp[1] += ftdict[(cidx,i)][3]
            tmp[2] += ftdict[(cidx,i)][4]
    sa = np.sum(tmp[0] * surfdict[cidx])
    tmp.append(sa)
    return tmp[1:]



    
def commRun(commlist, ftdict, surfdict, state = [160, 160]):
    """
    Input Commlist

    """
    print(len(commlist))
    print(len(surfdict))
    #print commlist
    nsnap = len(commlist)
    stat1 = {i:[] for i in range(nsnap)}
    stat2 = {i:[] for i in range(nsnap)}

    for cidx in commlist:
        stat1[cidx].append(len(commlist[cidx]))
        stat2[cidx].append(len(commlist[cidx]))
        out = processCommSnap(commlist[cidx], ftdict, surfdict, cidx)
        stat1[cidx].append(out)
        out = processCommSnapAggregate(commlist[cidx], ftdict, surfdict, cidx)
        stat2[cidx].append([out])
    
    #print stat1
    #print stat2

    cst = np.cumsum(state)
    step = 160
    statsum = []
    for idx, stat in enumerate([stat1, stat2]):
        out = []
        for i in range(nsnap):
            if stat[i][0] == 0:
                tmp = [i, 0, 0, 0, 0]
                out.append(tmp)
                continue
            tmp1 = [i, stat[i][0]]
            for j in stat[i][1]:
                tmp = tmp1 +  j
                out.append(tmp)
        head = "idx npkt space score sasa\n"
        #for i in out: print i
        with open('stat' + str(idx+1) + '.dat', 'w') as f:
            f.write(head)
            np.savetxt(f, out, fmt='%d')

        out = np.array(out)
        _, index =  np.unique(out[:,0], return_index=True)
        out = out[index,:]
        #print out[:,1]
        for i in cst:
            idx = np.where(out[(i-step):i,1] > 0)[0]
            idx =  idx + i - step
            statsum.append( [1.0 * len(idx) / step * 100.,  
                             np.average(out[idx, 2]), np.std(out[idx,2]),
                             np.average(out[idx, 3]), np.std(out[idx,3]),
                             np.average(out[idx, 4]), np.std(out[idx,4])])
            
    
    #print statsum
    with open('statsum.dat', 'w') as f:
        f.write('pop space spacesd score scoresd sasa sasasd\n')
        np.savetxt(f, statsum, fmt='%d')


def assignLigandOcc(lig, fn = 'cmdpocketByState_sum.txt'):
    """
    Given a binder and a list of cdpkt,
    For each cdpkt, for each snapshot, check if the binder is overlap with the binder.
    Criteria
    If any binder heavy atom is within 1.6 of pocket vertices_coords center.
    """
    from retired.AlphaSpace import pdbinfo
    from scipy import spatial

    ligcrd = []
    with open(lig, 'r') as f:
        for line in f:
            if pdbinfo.isAtom(line) == 1 and pdbinfo.isHydrogen(line) == 0:
                ligcrd.append(pdbinfo.coord(line))

    print(ligcrd)
    pdbidx = []
    with open(fn, 'r') as f:
        for line in f:
            if line[0] != 'P':
                if float(line.split(',')[-3]) > 0:
                    pdbidx.append(int(line.split(',')[0]))
    print(pdbidx)
    pdblist = []
    for i in pdbidx:
        pdblist.append(str(i).zfill(3)+'.pdb')
    print(pdblist)
    tree = spatial.cKDTree(ligcrd)

    df = np.zeros((len(pdblist), 160))
    for idx, pdb in enumerate(pdblist):
        dict = {i:[] for i in range(1,161,1)}
        with open(pdb, 'r') as f:
            for line in f:
                if pdbinfo.resn(line) == 'VNT' and pdbinfo.atmn(line) in  [' C  ',' O  '] :
                    dict[int(pdbinfo.resi(line))].append(pdbinfo.coord(line))

        print(dict)
        for i in range(1,161,1):
            if dict[i] != []:
                for j in dict[i]:
                    if tree.query_ball_point(j, 1.6) != []:
                        df[idx,i-1]=1
                        break
    
    for idx, pdb in enumerate(pdblist):
        print(df[idx])
        with open('tmp', 'w') as f:
            f.write(',ligocc\n')
            np.savetxt(f, df[idx], fmt=',%2d')
        os.system( 'paste ' + pdb[:-3] + "txt tmp > tmp2" )
        os.system( 'mv tmp2 ' + pdb[:-3] + "pst.txt" )
    os.system('rm tmp')

    avg = np.average(df,axis=1)

    nd = []
    idx = 0
    with open(fn, 'r') as f:
        for line in f:
            if line[0] != 'P':
                if float(line.split(',')[-3]) > 0:
                    nd.append(avg[idx])
                    idx += 1
                else:
                    nd.append(0.)
    print(nd)
    print(len(nd))
    with open('tmp', 'w') as f:
        f.write(',ligocc\n')
        print(nd)
        np.savetxt(f, np.array(nd), fmt=',%5.2f')
    os.system( 'paste ' + fn + " tmp > tmp2" )
    os.system( 'mv tmp2 ' + fn[:-3] + "pst.txt" )
    os.system('rm tmp')



