# -*- coding: utf-8 -*-

import numpy as np
import numpy
import scipy
from scipy import stats
from scipy.special import cbrt
import math
from operator import itemgetter
from collections import defaultdict
from scipy.spatial.distance import cdist
#from get_geom import _GetGeometricalCenter, GetAtomDistance, get_coordinates_from_pdb, _get_types_dict, get_type_coords

def ctd_alpha(Coordinates):
    ctd = np.array([np.mean(Coordinates,axis=0)])
    ctd_distance = cdist(ctd,Coordinates)[0]
    ctd_alpha={}
    ctd_alpha['ctd_1']=round(np.mean(ctd_distance),3)
    ctd_alpha['ctd_2']=round(np.std(ctd_distance),3)
    ctd_alpha['ctd_3']=round(cbrt(stats.skew(ctd_distance)),3)
    return ctd_alpha

def cst_alpha(Coordinates):
    ctd = np.array([np.mean(Coordinates,axis=0)])
    ctd_distance = cdist(ctd,Coordinates)[0]

    cst_index = np.argmin(ctd_distance)
    cst = np.array([Coordinates[cst_index]])
    cst_distance = cdist(cst,Coordinates)[0]
    
    cst_alpha={}
    cst_alpha['cst_1']=round(np.mean(cst_distance),3)
    cst_alpha['cst_2']=round(np.std(cst_distance),3)
    cst_alpha['cst_3']=round(cbrt(stats.skew(cst_distance)),3)
    return cst_alpha

def fct_alpha(Coordinates):
    ctd = np.array([np.mean(Coordinates,axis=0)])
    ctd_distance = cdist(ctd,Coordinates)[0]

    fct_index = np.argmax(ctd_distance)
    fct = np.array([Coordinates[fct_index]])
    fct_distance = cdist(fct,Coordinates)[0]
    
    fct_alpha={}
    fct_alpha['fct_1']=round(np.mean(fct_distance),3)
    fct_alpha['fct_2']=round(np.std(fct_distance),3)
    fct_alpha['fct_3']=round(cbrt(stats.skew(fct_distance)),3)
    return fct_alpha


def ftf_alpha(Coordinates):
    ctd = np.array([np.mean(Coordinates,axis=0)])
    ctd_distance = cdist(ctd,Coordinates)[0]

    fct_index = np.argmax(ctd_distance)
    fct = np.array([Coordinates[fct_index]])
    fct_distance = cdist(fct,Coordinates)[0]
    
    ftf_index = np.argmax(fct_distance)
    ftf = np.array([Coordinates[ftf_index]])
    ftf_distance = cdist(ftf,Coordinates)[0]
    
    ftf_alpha={}
    ftf_alpha['ftf_1']=round(np.mean(ftf_distance),3)
    ftf_alpha['ftf_2']=round(np.std(ftf_distance),3)
    ftf_alpha['ftf_3']=round(cbrt(stats.skew(ftf_distance)),3)
    return ftf_alpha

def Get_USR_alpha_beta(Coordinates):

    result={}
    result.update(ctd_alpha(Coordinates))
    result.update(cst_alpha(Coordinates))
    result.update(fct_alpha(Coordinates))
    result.update(ftf_alpha(Coordinates))

    return result

if __name__=="__main__":

    print(None)
