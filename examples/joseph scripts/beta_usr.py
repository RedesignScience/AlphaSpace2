# -*- coding: utf-8 -*-

import numpy as np
from scipy import stats
from scipy.special import cbrt
import math
from operator import itemgetter
from collections import defaultdict
from get_geom import _GetGeometricalCenter, GetAtomDistance, get_coordinates_from_pdb, get_types_dict, get_type_coords

def ctd_alpha(Coordinates):
    ctd=_GetGeometricalCenter(Coordinates)
    ctd_distance=[]
    for j in Coordinates:
        ctd_distance.append(GetAtomDistance(ctd,[float(j[2]),float(j[3]),float(j[4])]))
    ctd_alpha={}
    ctd_alpha['ctd_1']=round(np.mean(ctd_distance),3)
    ctd_alpha['ctd_2']=round(np.std(ctd_distance),3)
    ctd_alpha['ctd_3']=round(cbrt(stats.skew(ctd_distance)),3)
    return ctd_alpha

def cst_alpha(Coordinates):
    ctd=_GetGeometricalCenter(Coordinates)
    ctd_distance=[]
    for j in Coordinates:
        ctd_distance.append(GetAtomDistance(ctd,[float(j[2]),float(j[3]),float(j[4])]))
    cst_index=min(enumerate(ctd_distance), key=itemgetter(1))[0]
    cst=[float(Coordinates[cst_index][2]),float(Coordinates[cst_index][3]),float(Coordinates[cst_index][4])]
    cst_distance=[]
    for j in Coordinates:
        cst_distance.append(GetAtomDistance(cst,[float(j[2]),float(j[3]),float(j[4])]))
    cst_alpha={}
    cst_alpha['cst_1']=round(np.mean(cst_distance),3)
    cst_alpha['cst_2']=round(np.std(cst_distance),3)
    cst_alpha['cst_3']=round(cbrt(stats.skew(cst_distance)),3)
    return cst_alpha

def fct_alpha(Coordinates):
    ctd=_GetGeometricalCenter(Coordinates)
    ctd_distance=[]
    for j in Coordinates:
        ctd_distance.append(GetAtomDistance(ctd,[float(j[2]),float(j[3]),float(j[4])]))
    fct_index=max(enumerate(ctd_distance), key=itemgetter(1))[0]
    fct=[float(Coordinates[fct_index][2]),float(Coordinates[fct_index][3]),float(Coordinates[fct_index][4])]
    fct_distance=[]
    for j in Coordinates:
        fct_distance.append(GetAtomDistance(fct,[float(j[2]),float(j[3]),float(j[4])]))
    fct_alpha={}
    fct_alpha['fct_1']=round(np.mean(fct_distance),3)
    fct_alpha['fct_2']=round(np.std(fct_distance),3)
    fct_alpha['fct_3']=round(cbrt(stats.skew(fct_distance)),3)
    return fct_alpha


def ftf_alpha(Coordinates):
    ctd=_GetGeometricalCenter(Coordinates)
    ctd_distance=[]
    for j in Coordinates:
        ctd_distance.append(GetAtomDistance(ctd,[float(j[2]),float(j[3]),float(j[4])]))
    fct_index=max(enumerate(ctd_distance), key=itemgetter(1))[0]
    fct=[float(Coordinates[fct_index][2]),float(Coordinates[fct_index][3]),float(Coordinates[fct_index][4])]
    fct_distance=[]
    for j in Coordinates:
        fct_distance.append(GetAtomDistance(fct,[float(j[2]),float(j[3]),float(j[4])]))
    ftf_index=max(enumerate(fct_distance), key=itemgetter(1))[0]
    ftf=[float(Coordinates[ftf_index][2]),float(Coordinates[ftf_index][3]),float(Coordinates[ftf_index][4])]
    ftf_distance=[]
    for j in Coordinates:
        ftf_distance.append(GetAtomDistance(ftf,[float(j[2]),float(j[3]),float(j[4])]))
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

    print(res)
