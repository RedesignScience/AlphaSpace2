
# standard modules
import os
import pickle
import time

import numpy as np
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as scidist
from scipy.sparse import dok_matrix

from AS_fcn import calc_psimsa, calc_simsa, combine_prot_lig, genCommunityPocket, runAlphaSpace
# out module
from DTM_old.md_classes import MDpocket, RFCMDpocket, RFMDpocket, SSCluster
from retired.AlphaSpace import readPDB


def genMDpocket(pdblist, 
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
    ftdict = {}
    surfdict = {}
    snapidx = []
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
        os.system("tar -cvf pdb" + str(i) + ".tar --remove-files pdb" + str(i))
        os.chdir('..')
        
        ssp = ss.sub_pockets
        
        npocket_list.append(len(ssp))    
        x = fill_pocket(natom, ssp)
        surfdict[(i)] = pdb.asa_prot
        
        # vol_simp
        # p.vol_simp(ss.simplices, ss.prot_pdb)
        # vol_exc
        # p.vol_exc(ss.simplices,ss.prot_pdb,ss.asa_radius2)
        
        # for community
        pktComm = genCommunityPocket(pdb, ss)
        commlist[i] = pktComm

        for j,p in enumerate(ssp):
            snapidx.append((i,j))
            p.set_vol()
            p.set_vresid()
            ftdict[(i,j)] = [x[j], list(p.atm_cntr()), \
                list(p.vert_cntr()), p.space, p.score, p.dict_resid_vol]
        
        if lig_resid == []:
            for i in ss.resid:
                if i not in lig_resid:
                    lig_resid.append(i)
        
        if istime: 
            print(i, time.time() - t)
            t = time.time()
    os.chdir(olddir)
    
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
        clust_dist = 0.75 * max(zmat[:,2])
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
        pickle.dump(commlist, open('commlist.pkl', 'w'))
    
    ## surface state clustering
    #sscluster = getSSC(mdpocketlist)
    
    ## PPI modulation 
    #ppimod(mdpocketlist, lig_resid)
    
    return mdpocketlist


def fill_pocket(natom, pocket):
    """Convert pocket into binary atom array
    
    :param natom: Number of atoms
    :param pocket: Instance of pockets
    :returns: ndarray([len(pocekt), natom])
    :rtype: ndarray
    
    .. note:: If atom j in pocket i, assign (i,j) = 1.
    
    """
    
    npocket = len(pocket)
    pocket_atom_array = np.zeros([npocket, natom])
    for i in range(npocket):
        for j in pocket[i].atoms:
            pocket_atom_array[i,j] = 1.
    return pocket_atom_array



def genRandom(numsnap, step):
    """Random generated snapshot index
    """
    import random
    
    index = [i for i in range(numsnap)]
    sieve_index = []

    for i in range(int(numsnap/step)):
        k = random.choice(index)
        sieve_index.append(k)
        index.pop(k)
    sieve_index = sorted(sieve_index)
    return sieve_index


def genEven(numsnap, step):
    """Evenly spaced snapshot index
    """
    sieve_index = []
    
    for i in range(numsnap):
        if i % step == 0:
            sieve_index.append(i)
    return sieve_index


def genAvgArray(pocket_atom_array, cluster, nmdpocket):
    """Center of pocket by averaging all pockets in cluster
    
    :param pocket_atom_array: Pocket atom array
    :param cluster: cluster index
    :param nmdpocket: number of clusters
    :returns: average array for each cluster

    .. warning:: Not Tested
    
    """

    mdpocket_avg_array = [[] for i in range(nmdpocket)]
    
    ## count number pocket for each cluster for future use
    pocket_count = []
    
    for i in range(nmdpocket):
        pocket_count.append(cluster.count(i+1))
    
    for i,j in zip(pocket_atom_array, cluster):
        if mdpocket_avg_array[j-1] == []:
            mdpocket_avg_array[j-1] = i
        else:
            mdpocket_avg_array[j-1] += i
            
    for i in range(nmdpocket):
        mdpocket_avg_array[i] = mdpocket_avg_array[i] / pocket_count[i]
        
    return mdpocket_avg_array


def groupAvg(cluster, value, isSAwt = True):
    """Calculate the group average of value for each cluster 
    """
    
    value = np.array(value)
    cluster = np.array(cluster)
    avg = []
    avgidx = [] # this stores the index of the closest snapshot in cluster 

    for i in range(max(cluster)):
        tmpidx = np.where(cluster == (i+1))[0]
        avg.append(np.average(value[tmpidx]))
        if isSAwt:
            avgidx.append(tmpidx[np.argmax(value[tmpidx])])
        else:
            avgidx.append(tmpidx[np.argmin(value[tmpidx])])
    return np.round(np.array(avg),2), np.array(avgidx)


def getSSC(mdpocketlist, isCut = False, cutoff = 0.05,  
           ssm = 'volume',
           isSAwt = False, metric = 'euclidean',
           sscmethod = 'distance', sscdist = 0.5,
           isPlot = True):
    """Function to do surface state cluster
    """ 
    
    # number of snapshots from MD excluding snapshot of scan
    nsnap = mdpocketlist[0].nsnap
    # number of snapshots for scanning
    numscan = mdpocketlist[0].numscan
    # get number of dpkt
    ndtm = len(mdpocketlist)
    
    # initial matrix with (number of snapshot, number of mdpkt)
    dtmM = np.zeros([nsnap+numscan, ndtm])
    # assign mdpktM with values of pocket score
    for idx,p in enumerate(mdpocketlist):
        for sid,ss in enumerate(p.sslist):
            if len(ss.pocketidx) > 0:
                if ssm == 'score':
                    dtmM[sid][idx] = ss.score
                elif ssm == 'volume':
                    dtmM[sid][idx] = ss.volume

                else:           
                    dtmM[sid][idx] = 1               
    
    ## index of the dpkt for clustering
    ## if cut the outlier pocket based on the population of dpkt
    ## then cutoff will be applied and idx is the dpkt index used 
    ## for clustering
    idx = list(range(ndtm))
    if isCut:
        occplist = []
        for p in mdpocketlist:
            occplist.append(p.occprob)
        idx = np.where(np.array(occplist) > cutoff)[0]
        ndtmcut = len(list(idx))
    
    dtmM = dtmM[:,idx]
    # print dtmM
    
    if numscan > 0:
        dtmScan = dtmM[list(range(nsnap, nsnap+numscan)),]
        dtmM = dtmM[list(range(nsnap)),]
    
    # center distance
    cd =  np.average(np.sum(dtmM, axis=1))
    
    # pickle.dump(dtmM, open('dtmM.pkl', 'w'))
    
    ## pairwise distance metric
    ## this controls score weighted jaccard similarity
    if isSAwt:
        print('Metric is surface area weighted jaccard')
        data_dist = calc_psimsa(dtmM)
    elif metric == 'euclidean':
        print('Metric is euclidean')
        data_dist = scidist.pdist(dtmM)
    elif metric == 'jaccard':
        print('Metric is jaccard')
        data_dist = scidist.pdist(dtmM, "jaccard")
    else:
        print ("Error: wrong clustering method")
    
    # print data_dist.shape
    # print type(data_dist)
    # Z matrix
    zmat = hier.linkage(data_dist, method="average")

    ## plot heatmap with and without clustering
    if isPlot:
        plot_heatmap(data_dist, figname = "jaccard_cluster.png", isDist = True)
        plot_heatmap_nocluster(data_dist, figname = "jaccard_nocluster.png", isDist = True)
    
    ## cut the tree    
    if sscmethod == 'distance':
        clust_dist =  sscdist*max(zmat[:,2])
        #clust_dist = cd * 0.35
        print('DTM: Cutoff method distance, ', clust_dist)
        cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
    elif sscmethod == 'maxclust':
        clust_num = int(np.average(npocket_list)) 
        print('DTM: Cutoff method maxclust, ', clust_num)
        cluster = hier.fcluster(zmat,clust_num,criterion='maxclust')
    else:
        print("Error: wrong cutoff method")

    cluster = list(cluster)
    nssc = max(cluster)
    
    ## reindex the cluster id by the number of snapshot in each cluster
    clustercount = [cluster.count(i) for i in range(1, nssc+1,1)]  
    clustercountorder = list(np.argsort(clustercount)[::-1] + 1)
    clusterReorderdict = {j:i for i,j in zip(list(range(1,nssc+1,1)), clustercountorder)}
    cluster = [clusterReorderdict[i] for i in cluster]

    ## create SSCluster instance
    ssc = SSCluster(nssc) 
    ssc.cluster = cluster
    scoresum = np.sum(dtmM, axis=1)

    cluster = np.array(cluster)
    integrity = []
    
    data_dist_sqform = scidist.squareform(data_dist)

    for i in range(nssc):
        mat = np.outer(cluster == (i+1), cluster == (i+1))
        tmpidx = np.triu(mat, k=1)
        print(tmpidx)
        print(data_dist_sqform[tmpidx])
        print(type(data_dist_sqform[tmpidx]))
        if data_dist_sqform[tmpidx].shape[0] == 0:
            integrity.append(0.)        
        else:
            integrity.append(1.0 - np.average(data_dist_sqform[tmpidx]))
    

    # print integrity
    ssc.sscintegrity = integrity
    
    
    for p in mdpocketlist:
        p.getClustIntegrity(cluster)
        p.getClustScore(cluster)   
        #print p.clustintegrity
        ssc.sscclustintegrity.append(p.clustintegrity)
        ssc.sscclustscore.append(p.clustscore)
        ssc.sscclustvolume.append(p.clustvolume)
        ssc.sscclustpop.append(p.clustpop)
    cluster = list(cluster)
        
    for i in range(nssc):
        # get idx of snapshots is cluster i+1
        tmpidx = np.where(np.array(cluster) == (i+1))[0]
        
        # get the average score of cluster
        ssc.sscscoreavg.append(np.average(scoresum[tmpidx]).astype(int))
        
        # get the dtmM of cluster
        dtmMc = dtmM[tmpidx,]
        # get the average representation of the cluster
        tmp = np.average(dtmMc, axis=0)
        if isSAwt:
            # calculate the score weighted jaccard similarity between cluster and the average
            dtmMcdist = np.array([calc_simsa(dtmM[j,:], tmp) for j in tmpidx])
            # get the max similar index
            ssc.sscrep.append(tmpidx[np.argmax(dtmMcdist, axis = 0)])
        else:
            # calcuate the distance between cluster and the average
            dtmMcdist = np.sum((dtmMc-tmp)**2, axis=1)
            # get the min distance index
            ssc.sscrep.append(tmpidx[np.argmin(dtmMcdist, axis = 0)])

    if numscan > 0:
        ssc.scanlist(numscan)
        for i in range(numscan):
            if isSAwt:
                dtmMdist = np.array([ calc_simsa(dtmM[j,:], dtmScan[i,]) for j in range(nsnap)])
                idx = np.argmax(dtmMdist, axis=0)
            else:
                dtmMdist = np.sqrt(np.sum((dtmM-dtmScan[i,])**2, axis=1))
                idx = np.argmin(dtmMdist, axis=0)
            avg, avgidx = groupAvg(cluster, dtmMdist)
            ssc.assignScan(np.round(dtmMdist,2), avg, avgidx, cluster[idx], idx)

    ssc.writeSSC()
    ssc.writeSSCsum()
    ssc.writeIntegrity()
    
    # file for chimera 
    sscdict = {'snapshot': {i:{ j:[] for j in range(ndtm)} for i in range(nsnap+numscan)}}
    sscdict['cluster'] = {i:[] for i in range(nssc)}

    # all id are 0 indexed
    for idx,p in enumerate(mdpocketlist):
        for i,j in list(p.ftdict.keys()):
            #if i < nsnap:
            #print i, idx, j
            sscdict['snapshot'][i][idx].append(j)

    for idx,cidx in enumerate(cluster):
        sscdict['cluster'][cidx-1].append(idx)          
     
    sscdict['rep'] = ssc.sscrep
 
    sscdict['c'] = {i:{ j:[[],[]] for j in range(ndtm)} for i in range(nssc)}
    if isCut:
        sscdict['cut'] = ndtmcut
    else:
        sscdict['cut'] = -1
    
 
    for idx,p in enumerate(mdpocketlist):
        for sid, ss in enumerate(p.sslist):
            if ss.pocketidx == []:
                continue
            if sid < nsnap:
                sscdict['c'][cluster[sid]-1][idx][0].append(ss.atomcnt)
                sscdict['c'][cluster[sid]-1][idx][1].append(ss.vertcnt)    
 
    ssc.writeClusterPDB(sscdict['c'], ndtm)
    
    pickle.dump(sscdict, open('SSC_Chimera.pkl', 'w'))
    
    
def ppimod(mdpocketlist, lig_resid):
    """Function to get the modularity in ppi
    
    """
    nsnap = mdpocketlist[0].nsnap
    ndtm = len(mdpocketlist)
    nres = len(lig_resid)
    # print lig_resid

    dpkt_ss_res = {'ndtm' : ndtm, 'resid': lig_resid, 'nsnap': nsnap, 
                   'volmat' : { i : dok_matrix((nsnap, nres+1)) for i in range(ndtm)},
                   'occmat' : { i : dok_matrix((nsnap, nres+1)) for i in range(ndtm)},   # occ over total volume 
                   'occmat2' : { i : dok_matrix((nsnap, nres+1)) for i in range(ndtm)} } # occ over total occupied volume
                   #'volmat' : { i : 0 for i in range(ndtm)}}
    for idx,dp in enumerate(mdpocketlist):
        if not dp:
            continue
        for sid,ss in enumerate(dp.sslist):
            if len(ss.pocketidx) > 0:
                for r in ss.dict_resid_vol:
                    #print ss.dict_resid_vol
                    #print lig_resid.index(r)
                    dpkt_ss_res['volmat'][idx][sid, lig_resid.index(r)] = ss.dict_resid_vol[r]
                    dpkt_ss_res['occmat2'][idx][sid, lig_resid.index(r)] = ss.dict_resid_vol[r]
                    dpkt_ss_res['occmat'][idx][sid, lig_resid.index(r)] = ss.dict_resid_vol[r]/ss.volume
                dpkt_ss_res['volmat'][idx][sid,nres] = ss.volume - dpkt_ss_res['volmat'][idx][sid,::-1].sum()
                dpkt_ss_res['occmat'][idx][sid,nres] = dpkt_ss_res['volmat'][idx][sid,nres]/ss.volume
                dpkt_ss_res['occmat2'][idx][sid,::-1] = dpkt_ss_res['volmat'][idx][sid,::-1]/dpkt_ss_res['volmat'][idx][sid,::-1].sum()
                #print dpkt_ss_res['volmat'][idx][sid,nres], ss.volume, dpkt_ss_res['volmat'][idx][sid,nres]/ss.volume
                
    pickle.dump(dpkt_ss_res, open('dpkt_ss_res.pkl', 'w'))
    #print dpkt_ss_res
    

def Refine(mdpocketlist, srange=0):
    """Reclustering Dpkt with in srange by 
        
        D = [D1, D2, ...]
        Di = [[xS1, yS1, zS1], [xS2, yS2, zS2], ...]
        xSi = average([list of alpha sphere center x for Di in snapshot Si])
        DistM[i,j] is average pairwise distances between Dpkt i and j
        Average linkage clustering with cutoff 4

    :param srange:

    """

    # number of snapshots from MD excluding snapshot of scan
    nsnap = mdpocketlist[0].nsnap
    ndtm = len(mdpocketlist)  # get number of dpkt

    if srange == 0:
        srange = list(range(nsnap))
    
    dp = {i:i for i in range(ndtm)}

    print(srange)

    dtmM = [[] for i in range(ndtm)]
    for idx,p in enumerate(mdpocketlist):
        for sidx in srange:
            ss = p.sslist[sidx]
            if len(ss.pocketidx) > 0:
                if len(ss.vertcnt) > 1:
                    dtmM[idx].append(list(np.average(ss.vertcnt, axis=0)))
                dtmM[idx].append(ss.vertcnt[0])

    for i,j in enumerate(dtmM):
        if j == []:
            dp[i] = -1
    print(dp)

    distM = []
    idxlist = [ (i,j) for i in range(ndtm-1) for j in range(i+1, ndtm)]
    for i,j in idxlist:
        if dtmM[i] != [] and dtmM[j] != []: 
            d =  np.average(scidist.cdist(dtmM[i], dtmM[j]))
            distM.append(d)
            if d < 4:
                print(d, (i,j))

    zmat = hier.linkage(distM, method="average")
    rfcluster = hier.fcluster(zmat, 5, criterion='distance')
    print(rfcluster)

    k = 0
    for i in range(len(rfcluster)):
        while (dp[k] == -1):
            k = k + 1
        dp[k] = rfcluster[i]
        k = k + 1
    print(dp)

    rfcn = np.array([dp[i] for i in range(len(dp))])
    nrfcn = np.copy(rfcn)

    for i in range(1, max(rfcn)+1):
        loc = np.where(rfcn == i)
        nrfcn[loc] = loc[0][0]

    for i,j in enumerate(rfcn):
        print(i, j, nrfcn[i])
    print(nrfcn)

    return nrfcn
    

def StateSplitRefine(mdpocketlist, mdpocketidx, ftdict, surfdict, tmppdb, state = 0, isRefine = True, isCMD = False):
    """
    Function to split the DTM by State w/o Refine step
    based on the State index
    
    """
    # number of snapshots from MD excluding snapshot of scan
    nsnap = mdpocketlist[0].nsnap

    # get number of dpkt
    ndtm = len(mdpocketlist)

    if state == 0:
        srange = [list(range(nsnap))]
    else:
        state = np.cumsum(state)
        srange = [list(range(state[0]))]
        for i,j in enumerate(state[1:]):
            srange.append(list(range(state[i], state[i+1])))

    #print srange
    #print ns
    #print len(mdpocketidx)
    #print mdpocketidx

    ns = len(srange)
    rfdpktidx = [[[] for i in range(ndtm) ] for i in range(ns)]

    for i in range(ns):
        if isRefine:
            rfcluster = Refine(mdpocketlist, srange[i])
        else:
            rfcluster = list(range(ndtm))
        print(len(rfcluster))
        print(rfcluster)
        print(ndtm)
        for j, idxs in enumerate(mdpocketidx):
            for idx in idxs:
                if idx[0] in srange[i]:
                    rfdpktidx[i][rfcluster[j]].append(idx)

    #l = 0
    #for i,j,k in zip(rfcluster, rfdpktidx[0], mdpocketidx):
    #    print l, i,len(k), len(j)
    #    l += 1


    print(rfdpktidx)
    ## to do 
    ## new pocket assignment
    rfdpktlist = [[] for i in range(ns)]

    for i,idxs in enumerate(rfdpktidx):
        for idx in idxs:
            if idx == []:
                rfdpktlist[i].append(None)
            else:
                tmp = dict((k, ftdict[k]) for k in idx)
                if isCMD:
                    rfdpktlist[i].append(RFCMDpocket(tmp, idx, srange[i], surfdict))
                else:
                    rfdpktlist[i].append(RFMDpocket(tmp, idx, srange[i], surfdict))

    print(rfdpktlist)

    pickle.dump(rfdpktlist, open('rfdpktlist.pkl', 'w'))
    if isCMD:
        write_cmdpocketByState(rfdpktlist, tmppdb, state)
    else:
        write_mdpocketByState(rfdpktlist, tmppdb, state)


def write_mdpocket(mdpocketlist, pdb, nsnap):
    """Write the result from mdpocket in folder md_pockets
    
    :param mdpocketlist: list of pockets.
    :type mdpocketlist: str.
    :param pdb: original pdb
    :type pdb: str.
    :param nsnap: number of snapshot
    :type nsnap: int.
    :returns: folder of files
    
        - For each pocket, #.pdb have 
        
        -- atoms, center of atoms and center of alpha sphere
        
        - For each pocket, #.txt have
        
        -- frame #, # of pocket, total surf, total volume 
    
        - mdpocket_sum.txt have 
        
        -- pocket #, average p#, surf avg, surf std, volume avg, surf std
    
    .. todo:: nothing here
    """
    
    # Check folder md_pockets or create md_pockets
    if os.path.isdir("md_pockets"):
        os.system("rm md_pockets/*")
    else:
        os.system("mkdir md_pockets")

    # create file mdpocket_sum.txt
    mdsum = open("md_pockets/mdpocket_sum.txt","w")
    mdsum.write("""#PocketID  Pkt_avg  Pkt_std Surf_avg Surf_std  Vol_avg  Vol_std Score_avg Score_std Occupancy Integrity\n""")
    
    for idx,p in enumerate(mdpocketlist):
        ## write PDB for each pocket
        f = open("md_pockets/" + str(idx+1).zfill(3) + ".pdb", "w")
        for i,j in zip(p.atompop, pdb.prot):
            if i > 0:
                (j) = (j[0:60] + "%6.2f" + j[66:] + "\n") %(i)
                f.write(j)
        
        for sid, ss in enumerate(p.sslist):
            if ss.pocketidx == []:
                continue
            for crd in ss.atomcnt:
                f.write("""ATOM      1  C   CNT  %4d    """ %(sid+1))
                for i in crd:
                    f.write("%8.3f" %(i))
                f.write("\n")  
            for crd in ss.vertcnt:
                f.write("""ATOM      1  C   VNT  %4d    """ %(sid+1))
                for i in crd:
                    f.write("%8.3f" %(i))
                f.write("\n")  
        f.close()
    
        ## write Feature for each pocket
        g = open("md_pockets/" + str(idx+1).zfill(3) + ".txt", "w")
        g.write("""#Frame #pocket   SurfArea    Volume  Score\n""")
        
        for i,ss in enumerate(p.sslist):
            g.write("""%6d %6d %8.2f %8.2f %8.2f \n"""%(i+1, len(ss.pocketidx), ss.surf, ss.volume, ss.score))
        g.close()

        
        p.getIntegrity()

        mdsum.write("""%8d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n"""%(idx+1,p.avgstd_pktnum[0], \
                    p.avgstd_pktnum[1], p.avgstd_surf[0], p.avgstd_surf[1], p.avgstd_volume[0],\
                    p.avgstd_volume[1], p.avgstd_score[0],p.avgstd_score[1], p.occprob, p.integrity))
    
    mdsum.close()
    

def write_mdpocketByState(rfdpktlist, pdb, nsnap):
    """Write the result from mdpocket in folder md_pockets
    
    :param mdpocketlist: list of pockets.
    :type mdpocketlist: str.
    :param pdb: original pdb
    :type pdb: str.
    :param nsnap: number of snapshot
    :type nsnap: int.
    :returns: folder of files
    
        - For each pocket, #.pdb have 
        
        -- atoms, center of atoms and center of alpha sphere
        
        - For each pocket, #.txt have
        
        -- frame #, # of pocket, total surf, total volume 
    
        - mdpocket_sum.txt have 
        
        -- pocket #, average p#, surf avg, surf std, volume avg, surf std
    
    .. todo:: nothing here
    """
    
    # Check folder md_pockets or create md_pockets
    if os.path.isdir("md_pocketsByState"):
        os.system("rm md_pocketsByState/*/*")
    else:
        os.system("mkdir md_pocketsByState")


    ns = len(rfdpktlist)

    for stid in range(ns):
        mdpocketlist = rfdpktlist[stid]
        ssnap = nsnap[stid]

        if os.path.isdir("md_pocketsByState/" + str(stid)):
            os.system("rm md_pocketsByState/" + str(stid) + "/*")
        else:
            os.system("mkdir md_pocketsByState/" + str(stid))


        # create file mdpocket_sum.txt
        mdsum = open("md_pocketsByState/" + str(stid) + "/mdpocketByState_sum.txt","w")
        mdsum.write("""#PocketID  Pkt_avg  Pkt_std Surf_avg Surf_std  Vol_avg  Vol_std Score_avg Score_std Occupancy Integrity Diffusion\n""")
        
        for idx,p in enumerate(mdpocketlist):
            if p == None:
                mdsum.write("""%8d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n"""
                    %(idx+1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
                continue

            ## write PDB for each pocket
            f = open("md_pocketsByState/" + str(stid) + "/" + str(idx+1).zfill(3) + ".pdb", "w")
            for i,j in zip(p.atompop, pdb.prot):
                if i > 0:
                    (j) = (j[0:60] + "%6.2f" + j[66:] + "\n") %(i)
                    f.write(j)
            
            cntlist = []
            vntlist = [] 
            print(idx)
            for sid, ss in enumerate(p.sslist):
                if ss.pocketidx == []:
                    continue
                tmpcnt = []
                tmpvnt = []

                for crd in ss.atomcnt:
                    f.write("""ATOM      1  C   CNT  %4d    """ %(sid+1))
                    tmpcnt.append(crd)
                    for i in crd:
                        f.write("%8.3f" %(i))
                    f.write("\n")          
                
                for crd in ss.vertcnt:
                    f.write("""ATOM      1  C   VNT  %4d    """ %(sid+1))
                    tmpvnt.append(crd)
                    for i in crd:
                        f.write("%8.3f" %(i))
                    f.write("\n") 
                #print tmpcnt
                #print tmpvnt
                cntlist.append(np.average(tmpcnt, axis=0))
                vntlist.append(np.average(tmpvnt, axis=0))
            #print cntlist
            #print vntlist

            rad = np.sqrt(np.sum(np.std(vntlist, axis=0)**2))
            print(np.std(vntlist, axis=0))
            print(np.sqrt(np.sum(np.std(vntlist, axis=0)**2)))

            f.write("""ATOM      1  N   CNT  %4d    """ %(0))
            for i in np.average(cntlist, axis = 0):
                f.write("%8.3f" %(i))
            f.write("\n") 

            f.write("""ATOM      1  N   VNT  %4d    """ %(0))
            for i in np.average(vntlist, axis = 0):
                f.write("%8.3f" %(i))
            f.write("\n") 
            
            f.close()
        

            ## write Feature for each pocket
            g = open("md_pocketsByState/" + str(stid) + "/" + str(idx+1).zfill(3) + ".txt", "w")
            g.write("""#Frame #pocket   SurfArea    Volume  Score\n""")
            
            for i,ss in enumerate(p.sslist):
                g.write("""%6d %6d %8.2f %8.2f %8.2f \n"""%(i+1, len(ss.pocketidx), ss.surf, ss.volume, ss.score))
            g.close()

            
            p.getIntegrity()

            mdsum.write("""%8d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n"""%(idx+1,p.avgstd_pktnum[0], \
                        p.avgstd_pktnum[1], p.avgstd_surf[0], p.avgstd_surf[1], p.avgstd_volume[0],\
                        p.avgstd_volume[1], p.avgstd_score[0],p.avgstd_score[1], p.occprob, p.integrity, rad))
        
        mdsum.close()
    



def write_cmdpocketByState(rfdpktlist, pdb, nsnap):
    """Write the result from cmdpocket in folder cmd_pockets_byState
    
    :param mdpocketlist: list of pockets.
    :type mdpocketlist: str.
    :param pdb: original pdb
    :type pdb: str.
    :param nsnap: number of snapshot
    :type nsnap: int.
    :returns: folder of files
    
        - For each pocket, #.pdb have 
        
        -- atoms, center of atoms and center of alpha sphere
        
        - For each pocket, #.txt have
        
        -- frame #, # of pocket, total surf, total volume 
    
        - mdpocket_sum.txt have 
        
        -- pocket #, average p#, surf avg, surf std, volume avg, surf std
    
    .. todo:: nothing here
    """
    
    # Check folder md_pockets or create md_pockets
    if os.path.isdir("cmd_pocketsByState"):
        os.system("rm cmd_pocketsByState/*/*")
    else:
        os.system("mkdir cmd_pocketsByState")


    ns = len(rfdpktlist)

    for stid in range(ns):
        mdpocketlist = rfdpktlist[stid]
        ssnap = nsnap[stid]

        if os.path.isdir("cmd_pocketsByState/" + str(stid)):
            os.system("rm cmd_pocketsByState/" + str(stid) + "/*")
        else:
            os.system("mkdir cmd_pocketsByState/" + str(stid))


        # create file mdpocket_sum.txt
        mdsum = open("cmd_pocketsByState/" + str(stid) + "/cmdpocketByState_sum.txt","w")
        mdsum.write("""PocketID,Pkt_avg,Pkt_std,sasa_avg,sasa_sd,space_avg,space_sdd,score_avg,score_sd,coresasa_avg,coresasa_sd,corespace_avg,corespace_sdd,corescore_avg,corescore_sd,stability,integrity,dispersion\n""")
        

        outfmt = "%d," + ",".join(["%.2f"]  * 17 ) + "\n"
        for idx,p in enumerate(mdpocketlist):
            if p == None:
                mdsum.write(outfmt %(idx+1,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
                continue

            ## write PDB for each pocket
            f = open("cmd_pocketsByState/" + str(stid) + "/" + str(idx+1).zfill(3) + ".pdb", "w")
            for i,j in zip(p.atompop, pdb.prot):
                if i > 0:
                    (j) = (j[0:60] + "%6.2f" + j[66:] + "\n") %(i)
                    f.write(j)
            
            cntlist = []
            vntlist = [] 
            print(idx)
            
            for sid, ss in enumerate(p.sslist):
                if ss.pocketidx == []:
                    continue
                tmpcnt = []
                tmpvnt = []

                for crd, spc in zip(ss.atomcnt, ss.spacelist):
                    if spc >= 100.:
                        f.write("""ATOM      1  C   CNT  %4d    """ %(sid+1))
                        tmpcnt.append(crd)
                    else:
                        f.write("""ATOM      1  O   CNT  %4d    """ %(sid+1))

                    for i in crd:
                        f.write("%8.3f" %(i))
                    f.write("\n")      

                f.write("""ATOM      1  S   CNT  %4d    """ %(sid+1))
                for i in np.average(tmpcnt, axis=0):
                        f.write("%8.3f" %(i))
                f.write("\n")   



                for crd, spc in zip(ss.vertcnt, ss.spacelist):
                    if spc >= 100.:
                        f.write("""ATOM      1  C   VNT  %4d    """ %(sid+1))
                        tmpvnt.append(crd)
                    else:
                        f.write("""ATOM      1  O   VNT  %4d    """ %(sid+1))

                    for i in crd:
                        f.write("%8.3f" %(i))
                    f.write("\n") 


                f.write("""ATOM      1  S   VNT  %4d    """ %(sid+1))
                for i in np.average(tmpvnt, axis=0):
                        f.write("%8.3f" %(i))
                f.write("\n")  
                #print tmpcnt
                #print tmpvnt
                cntlist.append(np.average(tmpcnt, axis=0))
                vntlist.append(np.average(tmpvnt, axis=0))




            #print cntlist
            #print vntlist

            rad = np.sqrt(np.sum(np.std(vntlist, axis=0)**2))
            print(np.std(vntlist, axis=0))
            print(np.sqrt(np.sum(np.std(vntlist, axis=0)**2)))

            f.write("""ATOM      1  N   CNT  %4d    """ %(0))
            for i in np.average(cntlist, axis = 0):
                f.write("%8.3f" %(i))
            f.write("\n") 

            f.write("""ATOM      1  N   VNT  %4d    """ %(0))
            for i in np.average(vntlist, axis = 0):
                f.write("%8.3f" %(i))
            f.write("\n") 
            
            f.close()
        
            ## write Feature for each pocket
            g = open("cmd_pocketsByState/" + str(stid) + "/" + str(idx+1).zfill(3) + ".txt", "w")
            g.write("""frame,pocket,sasa,space,score,coresasa,corespace,corescore\n""")
            
            for i,ss in enumerate(p.sslist):
                g.write("""%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n"""
                        %(i+1, 
                        len(ss.pocketidx),
                        ss.surf, ss.volume, ss.score, 
                        ss.coresurf, ss.corevolume, ss.corescore))
            g.close()

            p.getIntegrity()
            mdsum.write(outfmt 
                        %(idx+1,
                        p.avgstd_pktnum[0], p.avgstd_pktnum[1], 
                        p.avgstd_surf[0], p.avgstd_surf[1], 
                        p.avgstd_volume[0], p.avgstd_volume[1], 
                        p.avgstd_score[0],p.avgstd_score[1], 
                        p.avgstd_coresurf[0], p.avgstd_coresurf[1], 
                        p.avgstd_corevolume[0], p.avgstd_corevolume[1], 
                        p.avgstd_corescore[0],p.avgstd_corescore[1],                     
                        p.occprob, p.integrity, rad))
        
        mdsum.close()
    




def plot_heatmap(pocket_atom_array,  metric="euclidean",
                 figname = "md_clust_heatmap.png", isDist = False):
    """Plot Heatmap of clustering
    
    This function is useful not only for this project and it is helpful to put 
    into other modules
    
    :param pocket_atom_array: the coord or the pairwised matrix in a list 
    :param metric: metric to calculated the distance, "euclidean" or "jaccard"
    :param figname: figure name
    :returns: figure 
    
    .. todo:: Should be put into other modules or plot modules
    
    """
    
    import scipy.spatial.distance as scidist 
    import scipy.cluster.hierarchy as hier
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.font_manager as font_manager
    
    # Font of figure
    path = '/Library/Frameworks/EPD64.framework/Versions/Current/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/Helvetica.ttf'
    try:
        prop = font_manager.FontProperties(fname=path)
        mpl.rcParams['font.family'] = prop.get_name()
    except:
        pass
    
    # color
    mpl.rc('axes',edgecolor='k')
    colors = ['#1e90ff', '#ff69b4', '#32cd32', '#9370db', '#2f4f4f','#87ceeb',\
          '#ff7f00','#228b22','#b22222','#000000']
    hier.set_link_color_palette(colors)
    
    # create pdist matrix (stored as list and see numpy for details)
    if isDist:
        data_dist = pocket_atom_array
    elif metric == 'euclidean':
        data_dist = scidist.pdist(pocket_atom_array)
    elif metric == 'jaccard':
        data_dist = scidist.pdist(pocket_atom_array, 'jaccard')
    else:
        print('Error')

    # figure size
    fig = plt.figure(figsize=(10,10))
    
    # x ywidth height
    ax1 = fig.add_axes([0.05,0.1,0.2,0.6])
    Y = hier.linkage(data_dist, method="average")
    Z1 = hier.dendrogram(Y, orientation='right',leaf_font_size = 10) 
    # adding/removing the axes
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.27,0.72,0.6,0.2])
    Z2 = hier.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    #Compute and plot the heatmap
    axmatrix = fig.add_axes([0.27,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = scidist.squareform(data_dist)
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=cm.RdBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([0.89,0.1,0.03,0.6])
    plt.colorbar(im, cax=axcolor)    
    plt.savefig(figname, dpi = 300,bbox_inches='tight')


def plot_heatmap_nocluster(pocket_atom_array,  
                 figname = "md_noclust_heatmap.png", isDist = False):
    """Plot Heatmap of clustering
    
    This function is useful not only for this project and it is helpful to put 
    into other modules
    
    :param pocket_atom_array: the coord or the pairwised matrix in a list 
    :param metric: metric to calculated the distance, "euclidean" or "jaccard"
    :param figname: figure name
    :returns: figure 
    
    .. todo:: Should be put into other modules or plot modules
    
    """
    
    import scipy.spatial.distance as scidist 
    import scipy.cluster.hierarchy as hier
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.font_manager as font_manager
    
    # Font of figure
    path = '/Library/Frameworks/EPD64.framework/Versions/Current/lib/python2.7/site-packages/matplotlib/mpl-data/fonts/ttf/Helvetica.ttf'
    try:
        prop = font_manager.FontProperties(fname=path)
        mpl.rcParams['font.family'] = prop.get_name()
    except:
        pass
    
    # color
    mpl.rc('axes',edgecolor='k')
    colors = ['#1e90ff', '#ff69b4', '#32cd32', '#9370db', '#2f4f4f','#87ceeb',\
          '#ff7f00','#228b22','#b22222','#000000']
    hier.set_link_color_palette(colors)
    
    # create pdist matrix (stored as list and see numpy for details)
    if isDist:
        data_dist = pocket_atom_array
    elif metric == 'euclidean':
        data_dist = scidist.pdist(pocket_atom_array)
    elif metric == 'jaccard':
        data_dist = scidist.pdist(pocket_atom_array, 'jaccard')
    else:
        print('Error')

    # figure size
    fig = plt.figure(figsize=(10,10))
    axmatrix = fig.add_axes([0.27,0.1,0.6,0.6])
    D = scidist.squareform(data_dist)
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=cm.RdBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([0.89,0.1,0.03,0.6])
    plt.colorbar(im, cax=axcolor)    
    plt.savefig(figname, dpi = 300,bbox_inches='tight')




def test1():
    """Test1: MDM2 Test"""
    pdblist = ['/Users/chengwang/Dropbox/Pocket_share/md/1ycr.pdb.'\
            + str(i) for i in range(1,79,1)]

    genMDpocket(pdblist, minsphere = 1, metric = 'jaccard',\
                cmethod = 'maxclust', sieve = False,  sieveopt = "random", \
                step = 5, isplot = True, isdebug = False, istime= True)

#    print mdpocketlist[0].atompop


def test2():
    """Test2: KIX Test for Surface State Clustering"""
    pdbs1 = [['../test/KIX/2AGH/0.' + str(i) + '.pdb'] for i in range(1,21)]
    pdbs2 = [['../test/KIX/2LXS/0.' + str(i) + '.pdb'] for i in range(1,21)]
    pdbs = pdbs1 + pdbs2
     
    
    mdpocketlist = genMDpocket(pdbs, minsphere =1, metric='jaccard', \
                 cmethod= 'maxclust', sieve=False, isplot=True, isdebug=False)


def test3():
    """Test3: KIX Test for Surface State Clustering with scanning of extra snapshot"""
    pdbs1 = [['../test/KIX/2AGH/0.' + str(i) + '.pdb'] for i in range(1,21)]
    pdbs2 = [['../test/KIX/2LXS/0.' + str(i) + '.pdb'] for i in range(1,21)]
    pdbs = pdbs1 + pdbs2
    scanpdbs = [pdbs2[19]]
     
    
    mdpocketlist = genMDpocket(pdbs, minsphere =1, metric='jaccard', \
                 cmethod= 'maxclust', sieve=False, isplot=True, isdebug=False,
                 isScan=True, scanlist = scanpdbs)


if __name__== "__main__":
    #test1()
    #test2()
    test3()
