
"""
Class for Community MD pocket:

CDTM
CMDpocket
CMDSnapshot

"""

import numpy as np
import scipy.spatial.distance as scidist
import os,sys

class CDTM:
    """
    Class for CDTM with information for each snapshot
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
        


class CMDpocket:
    """
    Class for each CMDpocket
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
    
class CMDSnapshot:
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

            
