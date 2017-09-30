
"""
Class for MD pocket:
MDpocket
MDSnapshot

"""

import numpy as np
import scipy.spatial.distance as scidist
import os,sys

class DTM:
    """
    Class for DTM with information for each snapshot
    """
    def __init__(self, nsnap):
        self.nsnap = nsnap
        self.idxlist = list(range(nsnap))
        self.sstate = []
        self.dtmM = []
        self.ndtm = 0
       
    def getdtmM(self, mdpocketlist):
        """
        Create the dtm matrix with dim nsnap, ndtm
        """
        self.ndtm = len(mdpocketlist)
        self.dtmM = np.zeros([self.nsnap, ndtm])
        for idx,p in enumerate(mdpocketlist):
            for sid, ss in enumerate(p.sslist):
                if len(ss.pocketidx) > 0:
                    dtmM[sid][idx] = ss.score
                    #dtmM[sid][idx] = 1               
                
        return self.ndtm, self.dtmM
        

class SSCluster:
    """
    Surface State Cluster class
    """
    def __init__(self, nssc):
        self.nssc = nssc
        self.cluster = []
        self.sscscoreavg = [] # average total score for each cluster
        self.sscintegrity = [] # average weighted similarity of each cluster (nssc)
        self.sscclustintegrity = [] # average jaccard similiarity of each dpkt within each 
                                    # cluster (nssc, ndpkt)
        self.sscclustscore = []
        self.sscclustvolume = []
        self.sscclustpop = []
        
        self.sscrep = []
        self.numscan = 0
        self.scandist = []
        self.scancidx = []
        self.scansidx = []
        self.scanavgdist = [] # this is a list with average distance(similarity)
        self.scanavgidx = []
        
    def scanlist(self, numscan):
        self.numscan = numscan
        
    def assignScan(self, dist, avg, avgidx,cidx, sidx):
        self.scandist.append(dist)
        self.scancidx.append(cidx)
        self.scansidx.append(sidx)
        self.scanavgdist.append(avg)
        self.scanavgidx.append(avgidx)
        
    def writeSSC(self, fn = 'SSC.dat'):
        outputM = [list(range(1, 1+ len(self.cluster)))]
        outputM.append(list(self.cluster))
        tmp = np.zeros(len(self.cluster))
        tmp[self.sscrep] = 1 
        outputM.append(list(tmp.astype(int).astype(str)))
        header = "snapshot  surfstat  represent "
        if self.numscan > 0:
            for i in range(self.numscan):
                outputM.append(list(np.round(self.scandist[i], 2).astype(str)))
                outputM.append([self.scancidx[i] for j in range(len(self.cluster))])
                header = header + "scan%-6sscand%-5s" %(str(i+1), str(i+1))
                
        outputM = np.transpose(np.array(outputM))
        f = open(fn, "w")
        f.write(header + "\n")
        for i in outputM:
            tmp = ''
            for j in i:
                tmp = tmp + '%-10s'  %(j) 
            tmp = tmp + "\n"
            f.write(tmp)
            
        f.close()
      
    def writeSSCsum(self, fn = 'SSC_sum.dat'):  
        header = 'surfstat  NumSS     Percent   RepID     AvgScore  '
        outputM = [list(range(1, 1 + self.nssc,1))]
        clustersnapshots = [self.cluster.count(i) for i in range(1, 1+self.nssc,1)]
        clustersnapshotspercent = list(np.round(1.0 * np.array(clustersnapshots)/len(self.cluster),2))
        outputM.append(clustersnapshots)
        outputM.append(clustersnapshotspercent)
        outputM.append(self.sscrep)
        # represent id is 0 indexed
        outputM.append(self.sscscoreavg)
        if self.numscan > 0:
            for i in range(self.numscan):
                outputM.append(list(self.scanavgdist[i].astype(str)))
                outputM.append(list(self.scanavgidx[i].astype(str)))
                header = header + "scan%-6sscanI%-5s" %(str(i+1), str(i+1))
        
        outputM = np.transpose(np.array(outputM))

        f = open(fn, "w")
        f.write(header + "\n")
        for i in outputM:
            tmp = ''
            for j in i:
                tmp = tmp + '%-10s'  %(j) 
            tmp = tmp + "\n"
            f.write(tmp)
            
        f.close()
        
    def writeClusterPDBold(self, sscdict, ndtm):
        """Write every pocket atom center and vertice center"""
        os.system('mkdir ssc')
        for i in range(self.nssc):
            fn = 'ssc/ssc' + str(i) + '.pdb'
            f = open(fn, 'w')
            
            for j in range(ndtm):
                if sscdict[i][j][0] != []:
                    for crds in sscdict[i][j][0]:
                        for crd in crds:
                            f.write("""ATOM      1  C   CNT  %4d    """ %(j+1))
                            for c in crd:
                                f.write("%8.3f" %(c))
                            f.write("\n")                        
                    for crds in sscdict[i][j][1]:
                        for crd in crds:
                            f.write("""ATOM      1  C   VNT  %4d    """ %(j+1))
                            for c in crd:
                                f.write("%8.3f" %(c))
                            f.write("\n") 
     
            f.close()


    def writeClusterPDB(self, sscdict, ndtm):
        """Write every pocket atom center and vertice center
            if more than one for each dtm pocket, it will be 
            averaged    
        """
        os.system('mkdir ssc')
        for i in range(self.nssc):
            fn = 'ssc/ssc' + str(i) + '.pdb'
            f = open(fn, 'w')
            
            for j in range(ndtm):
                if sscdict[i][j][0] != []:
                    for crds in sscdict[i][j][0]:
                        crd = np.average(crds, axis=0)

                        f.write("""ATOM      1  C   CNT  %4d    """ %(j+1))
                        for c in crd:
                            f.write("%8.3f" %(c))
                        f.write("\n")      
                                          
                    for crds in sscdict[i][j][1]:
                        crd = np.average(crds, axis=0)

                        f.write("""ATOM      1  C   VNT  %4d    """ %(j+1))
                        for c in crd:
                            f.write("%8.3f" %(c))
                        f.write("\n") 

            f.close()


    def writeIntegrity(self):
        """Write the integrity for each dpkt within each cluster
        
        """

        np.savetxt('ssc/SSC_integrity.txt', np.array(self.sscclustintegrity), fmt='%4.2f')
        np.savetxt('ssc/SSC_score.txt', np.array(self.sscclustscore), fmt='%8.2f')
        np.savetxt('ssc/SSC_volume.txt', np.array(self.sscclustvolume), fmt='%8.2f')
        np.savetxt('ssc/SSC_pop.txt', np.array(self.sscclustpop), fmt='%4.2f')

  

class MDpocket:
    """
    Class for each MDpocket
    """
    
    def __init__(self, ftdict, snapidx, nsnap, surfdict, numscan=0):
       
        # input variable
        self.ftdict = ftdict
        self.snapidx = snapidx
        self.nsnap = nsnap
        self.surfdict = surfdict
        self.avgstd_volume = []
        self.avgstd_surf = []
        self.avgstd_score = []
        self.numscan = numscan # number of snapshot for scanning
        self.integrity = -1.
        self.clustintegrity = [] # this is for surface state cluster
        self.clustscore = []
        self.clustvolume = []
        self.clustpop = []
        
        # get number of atoms
        natom = len(self.surfdict[0])
        self.atompop = np.zeros(natom) # list or array of atom populations   
        self.sslist = [MDSnapshot(i) for i in range(self.nsnap + self.numscan)]
               
        for i,j in enumerate(self.snapidx):   
            tmp1, tmp2 = j
            self.sslist[tmp1].set_pocketidx(j)
            self.sslist[tmp1].set_atomarray(self.ftdict[j][0].astype(int))      
            self.sslist[tmp1].set_atomcnt(self.ftdict[j][1])
            self.sslist[tmp1].set_vertcnt(self.ftdict[j][2])  
            self.sslist[tmp1].set_volume(self.ftdict[j][3])
            self.sslist[tmp1].set_score(self.ftdict[j][4]) ## new dict term
            self.sslist[tmp1].set_dict_resid_vol(self.ftdict[j][5]) ## dict of {resid:occvol}
        
        # assign surface area
        for i in range(self.nsnap + self.numscan):
            tmparray = self.sslist[i].atomarray
            if tmparray != []:
                self.sslist[i].set_surf(sum(tmparray * self.surfdict[(i)]))
        
        # only consider the snapshots in md not the scanning ones
        volumelist = []
        scorelist = []
        surflist = []
        pktnumlist = []
        for i in range(nsnap):
            tmparray = self.sslist[i].atomarray
            if tmparray != []:
                self.atompop += tmparray
                volumelist.append(self.sslist[i].volume)
                scorelist.append(self.sslist[i].score)
                surflist.append(self.sslist[i].surf)
                pktnumlist.append(len(self.sslist[i].pocketidx))
                
        self.avgstd_volume = [np.average(volumelist), np.std(volumelist)]
        self.avgstd_score = [np.average(scorelist), np.std(scorelist)]
        self.avgstd_surf = [np.average(surflist), np.std(surflist)]
        self.avgstd_pktnum = [np.average(pktnumlist), np.std(pktnumlist)]
        self.atompop = self.atompop / self.nsnap
        self.occprob = float(len(volumelist)) / self.nsnap
        
   
    def getIntegrity(self):
        """
        Calculate the integrity of the dpkt
        """ 
        atmarray = []  
        for i in range(self.nsnap):
            tmparray = self.sslist[i].atomarray
            if tmparray != []:
                atmarray.append(tmparray)
        if len(atmarray) == 1:
            self.integrity = 0.0
        else:
            self.integrity = 1.0 - np.average(scidist.pdist(np.array(atmarray), 'jaccard'))
   
    def getClustIntegrity(self, cluster):
        """
        Calculate the integrity of the dpkt
        """ 
        nssc = max(cluster)
        atmarray = []  
        for i in range(nssc):
            idxs = np.where(cluster == (i+1))[0]
            if idxs.shape[0] == 0:
                self.clustintegrity.append(-1)
            else:
                atmarray = []
                for j in idxs:                
                    tmparray = self.sslist[j].atomarray
                    if tmparray != []:
                        atmarray.append(tmparray)
                if atmarray == []:
                    self.clustintegrity.append(-1)
                elif len(atmarray) == 1:
                    self.clustintegrity.append(0)
                else:
                    self.clustintegrity.append(1.0-np.average(scidist.pdist(np.array(atmarray), 'jaccard')))


    def getClustScore(self, cluster):
        """
        Calculate the integrity of the dpkt
        """ 
        nssc = max(cluster)

        for i in range(nssc):
            idxs = np.where(cluster == (i+1))[0]
            if idxs.shape[0] == 0:
                self.clustscore.append(-1)
                self.clustvolume.append(-1)
                self.clustpop.append(0)
            else:
                scorearray = []
                volumearray = []
                
                for j in idxs:                
                    tmparray = self.sslist[j].atomarray
                    if tmparray != []:
                        scorearray.append(self.sslist[j].score)
                        volumearray.append(self.sslist[j].volume)
                if scorearray == []:
                    self.clustscore.append(-1)
                    self.clustvolume.append(-1)
                    self.clustpop.append(0)
                else:
                    self.clustscore.append(np.average(scorearray))
                    self.clustvolume.append(np.average(volumearray))  
                    self.clustpop.append(1.0 * len(scorearray)/len(cluster))




class RFMDpocket:
    """
    Class for each Refined MDpocket
    """
    
    def __init__(self, ftdict, snapidx, srange, surfdict):
       
        # input variable
        self.ftdict = ftdict
        self.snapidx = snapidx
        self.srange = srange
        self.nsnap = len(self.srange)
        self.surfdict = surfdict
        # three propeties volume, surf, score
        # list with [avg, sd]
        self.avgstd_volume = []
        self.avgstd_surf = []
        self.avgstd_score = []
        self.integrity = -1.
        self.clustintegrity = [] # this is for surface state cluster
        self.clustscore = []
        self.clustvolume = []
        self.clustpop = []
        
        # get number of atoms
        natom = len(self.surfdict[0])
        self.atompop = np.zeros(natom) # list or array of atom populations   
        self.sslist = [MDSnapshot(i) for i in self.srange]

        self.id2sid = { j:i for i,j in enumerate(self.srange)}

        for i,j in enumerate(self.snapidx):   
            tmp1, tmp2 = j
            tmp1 = self.id2sid[tmp1]
            self.sslist[tmp1].set_pocketidx(j)
            self.sslist[tmp1].set_atomarray(self.ftdict[j][0].astype(int))      
            self.sslist[tmp1].set_atomcnt(self.ftdict[j][1])
            self.sslist[tmp1].set_vertcnt(self.ftdict[j][2])  
            self.sslist[tmp1].set_volume(self.ftdict[j][3])
            self.sslist[tmp1].set_score(self.ftdict[j][4]) ## new dict term
            self.sslist[tmp1].set_dict_resid_vol(self.ftdict[j][5]) ## dict of {resid:occvol}
        
        # assign surface area
        for i in self.srange:
            id = self.id2sid[i]
            tmparray = self.sslist[id].atomarray
            if tmparray != []:
                self.sslist[id].set_surf(sum(tmparray * self.surfdict[(i)]))
        
        # only consider the snapshots in md not the scanning ones
        volumelist = []
        scorelist = []
        surflist = []
        pktnumlist = []
        for i in self.srange:
            id = self.id2sid[i]
            tmparray = self.sslist[id].atomarray
            if tmparray != []:
                self.atompop += tmparray
                volumelist.append(self.sslist[id].volume)
                scorelist.append(self.sslist[id].score)
                surflist.append(self.sslist[id].surf)
                pktnumlist.append(len(self.sslist[id].pocketidx))
                
        self.avgstd_volume = [np.average(volumelist), np.std(volumelist)]
        self.avgstd_score = [np.average(scorelist), np.std(scorelist)]
        self.avgstd_surf = [np.average(surflist), np.std(surflist)]
        self.avgstd_pktnum = [np.average(pktnumlist), np.std(pktnumlist)]
        self.atompop = self.atompop / self.nsnap
        self.occprob = float(len(volumelist)) / self.nsnap
        
    def getIntegrity(self):
        """
        Calculate the integrity of the dpkt
        """ 
        atmarray = []  
        for i in range(self.nsnap):
            tmparray = self.sslist[i].atomarray
            if tmparray != []:
                atmarray.append(tmparray)
        if len(atmarray) == 1:
            self.integrity = 0.0
        else:
            self.integrity = 1.0 - np.average(scidist.pdist(np.array(atmarray), 'jaccard'))



class RFCMDpocket:
    """
    Class for each Refined MDpocket
    """
    
    def __init__(self, ftdict, snapidx, srange, surfdict):
       
        # input variable
        self.ftdict = ftdict
        self.snapidx = snapidx
        self.srange = srange
        self.nsnap = len(self.srange)
        self.surfdict = surfdict
        
        # different from RFMDpocket 
        # three propeties volume, surf, score
        # list with [avg, sd]
        self.avgstd_volume = []
        self.avgstd_surf = []
        self.avgstd_score = []
        self.avgstd_corevolume = []
        self.avgstd_coresurf = []
        self.avgstd_corescore = []

        # integrity is calculated based on core pocket 
        self.integrity = -1.
        
        # this is for surface state cluster
        self.clustintegrity = []
        self.clustscore = []
        self.clustvolume = []
        self.clustpop = []
        
        # get number of atoms
        natom = len(self.surfdict[0])
        self.atompop = np.zeros(natom) # list or array of atom populations   
        self.sslist = [CMDSnapshot(i) for i in self.srange]

        self.id2sid = { j:i for i,j in enumerate(self.srange)}

        for i,j in enumerate(self.snapidx):   
            tmp1, tmp2 = j
            tmp1 = self.id2sid[tmp1]
            self.sslist[tmp1].set_pocketidx(j)
            if ftdict[j][3] >= 100.:
                self.sslist[tmp1].set_corevolume(self.ftdict[j][3])
                self.sslist[tmp1].set_corescore(self.ftdict[j][4]) ## new dict term
                self.sslist[tmp1].set_coreatomarray(self.ftdict[j][0].astype(int))
            
            self.sslist[tmp1].set_atomarray(self.ftdict[j][0].astype(int))      
            self.sslist[tmp1].set_atomcnt(self.ftdict[j][1])
            self.sslist[tmp1].set_vertcnt(self.ftdict[j][2])  
            self.sslist[tmp1].set_volume(self.ftdict[j][3])
            self.sslist[tmp1].set_score(self.ftdict[j][4]) ## new dict term
            self.sslist[tmp1].set_dict_resid_vol(self.ftdict[j][5]) ## dict of {resid:occvol}
        
        # assign surface area
        for i in self.srange:
            id = self.id2sid[i]
            tmparray = self.sslist[id].atomarray
            if tmparray != []:
                self.sslist[id].set_surf(sum(tmparray * self.surfdict[(i)]))
            coretmparray = self.sslist[id].coreatomarray
            if coretmparray != []:
                self.sslist[id].set_coresurf(sum(coretmparray * self.surfdict[(i)]))
        
        # only consider the snapshots in md not the scanning ones
        volumelist = []
        scorelist = []
        surflist = []
        corevolumelist = []
        corescorelist = []
        coresurflist = []
        pktnumlist = []
        for i in self.srange:
            id = self.id2sid[i]
            tmparray = self.sslist[id].atomarray
            if tmparray != []:
                self.atompop += tmparray
                volumelist.append(self.sslist[id].volume)
                scorelist.append(self.sslist[id].score)
                surflist.append(self.sslist[id].surf)
                corevolumelist.append(self.sslist[id].corevolume)
                corescorelist.append(self.sslist[id].corescore)
                coresurflist.append(self.sslist[id].coresurf)
                pktnumlist.append(len(self.sslist[id].pocketidx))
                
        self.avgstd_volume = [np.average(volumelist), np.std(volumelist)]
        self.avgstd_score = [np.average(scorelist), np.std(scorelist)]
        self.avgstd_surf = [np.average(surflist), np.std(surflist)]
        self.avgstd_corevolume = [np.average(corevolumelist), np.std(corevolumelist)]
        self.avgstd_corescore = [np.average(corescorelist), np.std(corescorelist)]
        self.avgstd_coresurf = [np.average(coresurflist), np.std(coresurflist)]

        self.avgstd_pktnum = [np.average(pktnumlist), np.std(pktnumlist)]
        self.atompop = self.atompop / self.nsnap
        self.occprob = float(len(volumelist)) / self.nsnap
        
    def getIntegrity(self):
        """
        Calculate the integrity of the dpkt
        """ 
        atmarray = []  
        for i in range(self.nsnap):
            tmparray = self.sslist[i].coreatomarray
            if tmparray != []:
                atmarray.append(tmparray)
        if len(atmarray) == 1:
            self.integrity = 0.0
        else:
            self.integrity = 1.0 - np.average(scidist.pdist(np.array(atmarray), 'jaccard'))



class MDSnapshot:
    """
    Pocket information for each snapshot
    """
    
    def __init__(self, ssid):
        self.ssid = ssid
        self.atomcnt = []
        self.vertcnt = []
        self.pocketidx = []
        self.surf = 0.
        self.volume = 0.
        self.spacelist = [] # this is for the community pocket aux and core output
        self.atomarray = []
        self.score = 0.
        self.dict_resid_vol = {}
        
    def set_atomcnt(self, atomcnt):
        self.atomcnt.append(atomcnt)
        
    def set_vertcnt(self, vertcnt):
        self.vertcnt.append(vertcnt)

    def set_pocketidx(self, pocketidx):
        self.pocketidx.append(pocketidx)
        
    def set_surf(self, surf):
        # not additive
        self.surf = surf
    
    def set_volume(self, volume):
        # additive
        self.volume += volume
        self.spacelist.append(volume)
        
    def set_score(self, score):
        # additive
        self.score += score
        
    def set_atomarray(self, atomarray):
        if self.atomarray == []:
            self.atomarray = atomarray
        else:
            self.atomarray = self.atomarray | atomarray
    
    def set_dict_resid_vol(self, dict_resid_vol):
        if self.dict_resid_vol == {}:
            self.dict_resid_vol = dict_resid_vol
        else:
            for i in dict_resid_vol:
                if i not in self.dict_resid_vol:
                    self.dict_resid_vol[i] = dict_resid_vol[i]
                else:
                    self.dict_resid_vol[i] += dict_resid_vol[i]


class CMDSnapshot:
    """
    Pocket information for each snapshot
    """
    
    def __init__(self, ssid):
        self.ssid = ssid
        self.atomcnt = []
        self.vertcnt = []
        self.pocketidx = []
        
        # core + aux
        self.surf = 0.
        self.volume = 0.
        self.score = 0.
        self.atomarray = []

        # core pocket property
        self.coresurf = 0.
        self.corevolume = 0.
        self.corescore = 0.
        self.coreatomarray = []

        # this is for the community pocket aux and core output
        self.spacelist = [] 
        self.atomarray = []
        self.dict_resid_vol = {}
        
    def set_atomcnt(self, atomcnt):
        self.atomcnt.append(atomcnt)
        
    def set_vertcnt(self, vertcnt):
        self.vertcnt.append(vertcnt)

    def set_pocketidx(self, pocketidx):
        self.pocketidx.append(pocketidx)
        
    def set_surf(self, surf):
        # not additive
        self.surf = surf
    
    def set_volume(self, volume):
        # additive
        self.volume += volume
        self.spacelist.append(volume)
        
    def set_score(self, score):
        # additive
        self.score += score

    def set_coresurf(self, coresurf):
        # not additive
        self.coresurf = coresurf
    
    def set_corevolume(self, corevolume):
        # additive
        self.corevolume += corevolume
        
    def set_corescore(self, corescore):
        # additive
        self.corescore += corescore

    def set_atomarray(self, atomarray):
        #
        if self.atomarray == []:
            self.atomarray = atomarray
        else:
            self.atomarray = self.atomarray | atomarray

    def set_coreatomarray(self, coreatomarray):
        #
        if self.coreatomarray == []:
            self.coreatomarray = coreatomarray
        else:
            self.coreatomarray = self.coreatomarray | coreatomarray
    
    def set_dict_resid_vol(self, dict_resid_vol):
        #
        if self.dict_resid_vol == {}:
            self.dict_resid_vol = dict_resid_vol
        else:
            for i in dict_resid_vol:
                if i not in self.dict_resid_vol:
                    self.dict_resid_vol[i] = dict_resid_vol[i]
                else:
                    self.dict_resid_vol[i] += dict_resid_vol[i]         
