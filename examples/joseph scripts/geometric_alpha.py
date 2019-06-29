# -*- coding: utf-8 -*-

import glob
import itertools
import numpy as np
import numpy
import scipy
from scipy import stats
from scipy.special import cbrt
import math
from operator import itemgetter
from AtomProperty import GetRelativeAtomicProperty
from collections import defaultdict
from get_geom import _GetGeometricalCenter, GetAtomDistance, get_coordinates_from_pdb, get_types_dict, GetGementricalDistanceMatrix, get_type_coords

def _GetMassCenter(MassCoordinates):
    res1=0.0
    res2=0.0
    res3=0.0
    temp=[]
    for i in MassCoordinates:
        res1=res1+i[0]*i[1][0]
        res2=res2+i[0]*i[1][1]
        res3=res3+i[0]*i[1][2]
        temp.append(i[0])
    result=[res1/sum(temp),res2/sum(temp),res3/sum(temp)]
    return result


def Calculate3DWiener(DistanceMatrix):

    return round((scipy.sum(DistanceMatrix)/2.0)/1000,3)

def CalculatePetitjean3DIndex(DistanceMatrix):

    temp1=scipy.amax(DistanceMatrix,axis=0)

    temp2=round(max(temp1)/min(temp1)-1.0,3)
    if math.isnan(temp2):
        return 0
    else:
        return temp2

def CalculateGemetricalDiameter(DistanceMatrix):

    temp1=scipy.amax(DistanceMatrix,axis=0)

    return round(max(temp1),3)

def CalculateGravitational3D1(Coordinates,types_dict):
    coords=[]
    for j in Coordinates:
        coords.append([float(types_dict[j[0]][j[1]][3]),[float(j[2]),float(j[3]),float(j[4])]])

    nAT=len(coords)
    result=0.0
    for i in range(nAT-1):
        for j in range(i+1,nAT):
            dis=GetAtomDistance(coords[i][1],coords[j][1])

            result=result+coords[i][0]*coords[j][0]/scipy.power(dis,p=2)

    return round(float(result)/100,3)

def CalculateRadiusofGyration(Coordinates,types_dict):
    coords=[]
    mass=0.0
    for j in Coordinates:
        coords.append([float(types_dict[j[0]][j[1]][3]),[float(j[2]),float(j[3]),float(j[4])]])
        mass=mass+float(types_dict[j[0]][j[1]][3])
    
    nAT=len(coords)
    masscenter=_GetMassCenter(coords)
    result=0.0
    for i in range(nAT):
        dis=GetAtomDistance(coords[i][1],masscenter)
        result=result+coords[i][0]*scipy.power(dis,p=2)


    return round(scipy.sqrt(float(result/mass)),3)

def CalculateHarary3D(DistanceMatrix):
    nAT=len(DistanceMatrix)
    res=0.0
    for i in range(nAT-1):
        for j in range(i+1,nAT):
            if DistanceMatrix[i,j]==0:
                cds=0.0
            else:
                cds=1./DistanceMatrix[i,j]
            res=res+cds
    return round(res/10,3)


def CalculateAverageGeometricalDistanceDegree(DistanceMatrix):
    nAT=len(DistanceMatrix)

    res=sum(sum(DistanceMatrix))/nAT

    return round(res,3)

def CalculateAbsEigenvalueSumOnGeometricMatrix(DistanceMatrix):
    u,s,vt=scipy.linalg.svd(DistanceMatrix)

    return round((sum(abs(s)))/10,3)

def CalculateSPANR(Coordinates,types_dict):
    coords=[]
    for j in Coordinates:
        coords.append([float(types_dict[j[0]][j[1]][3]),[float(j[2]),float(j[3]),float(j[4])]])

    masscenter=_GetMassCenter(coords)

    res=[]
    for i in coords:
        res.append(GetAtomDistance(i[1],masscenter))

    return round(float(max(res)),3)

def CalculateAverageSPANR(Coordinates,types_dict):
    coords=[]
    for j in Coordinates:
        coords.append([float(types_dict[j[0]][j[1]][3]),[float(j[2]),float(j[3]),float(j[4])]])

    masscenter=_GetMassCenter(coords)
    nAT=len(coords)

    res=[]
    for i in coords:
        res.append(GetAtomDistance(i[1],masscenter))

    return round(math.pow(float(max(res))/nAT,0.5),3)

def GetInertiaMatrix(Coordinates,types_dict):
    coords=[]
    mass_center=_GetGeometricalCenter(Coordinates)
    for j in Coordinates:
        coords.append([float(types_dict[j[0]][j[1]][3]),list(np.array([float(j[2]),float(j[3]),float(j[4])])-np.array(mass_center))])

    nAT=len(coords)
    
    InertiaMatrix=scipy.zeros((3,3))
    res11=0.0
    res22=0.0
    res33=0.0
    res12=0.0
    res23=0.0
    res13=0.0
    for i in range(nAT):
        res11=res11+coords[i][0]*(math.pow(coords[i][1][1],2)+math.pow(coords[i][1][2],2))
        res22=res22+coords[i][0]*(math.pow(coords[i][1][0],2)+math.pow(coords[i][1][2],2))
        res33=res33+coords[i][0]*(math.pow(coords[i][1][0],2)+math.pow(coords[i][1][1],2))
        res12=res12+coords[i][0]*(coords[i][1][0]*coords[i][1][1])
        res13=res13+coords[i][0]*(coords[i][1][0]*coords[i][1][2])
        res23=res23+coords[i][0]*(coords[i][1][1]*coords[i][1][2])
    InertiaMatrix[0,0]=res11
    InertiaMatrix[1,1]=res22
    InertiaMatrix[2,2]=res33
    InertiaMatrix[0,1]=-res12
    InertiaMatrix[0,2]=-res13
    InertiaMatrix[1,2]=-res23
    InertiaMatrix[1,0]=-res12
    InertiaMatrix[2,0]=-res13
    InertiaMatrix[2,1]=-res23

    return InertiaMatrix


def CalculateMolecularEccentricity(Coordinates,types_dict):

    InertiaMatrix=GetInertiaMatrix(Coordinates,types_dict)
    u,s,v=scipy.linalg.svd(InertiaMatrix)

    res1=s[0]
    res3=s[2]

    res=math.pow(res1*res1-res3*res3,1./2)/res1
    if math.isnan(res):
        return 0
    else:
        return round(res,3)

def CalculatePrincipalMomentofInertia(Coordinates,types_dict):
    InertiaMatrix=GetInertiaMatrix(Coordinates,types_dict)
    u,s,v=scipy.linalg.svd(InertiaMatrix)

    res={}
    res['IA']=round(s[2],3)
    res['IB']=round(s[1],3)
    res['IC']=round(s[0],3)

    return res

def CalculateRatioPMI(Coordinates,types_dict):
    temp=CalculatePrincipalMomentofInertia(Coordinates,types_dict)
    res={}
    if temp['IB'] == 0 or temp['IC'] == 0:
        res['IA/B']=0
        res['IA/C']=0
        res['IB/C']=0
    else:
        res['IA/B']=round(temp['IA']/temp['IB'],3)
        res['IA/C']=round(temp['IA']/temp['IC'],3)
        res['IB/C']=round(temp['IB']/temp['IC'],3)
    return res

def GetGeometric_withtypes(Coordinates,types_dict,atom_type):
    if atom_type not in ('HP','PL','N','P','D','A','DA','AR','no_type','HX'):
        return 'invalid atom type'
    elif atom_type in ('HP','PL','N','P','D','A','DA','AR','HX'):
        types_coords=get_type_coords(Coordinates,types_dict,atom_type)
        if types_coords:
            DistanceMatrix=GetGementricalDistanceMatrix(types_coords)
            res={}
            res['W3D_'+atom_type]=Calculate3DWiener(DistanceMatrix)
            res['Petitj3D_'+atom_type]=CalculatePetitjean3DIndex(DistanceMatrix)
            res['GeDi_'+atom_type]=CalculateGemetricalDiameter(DistanceMatrix)
            res['grav_'+atom_type]=CalculateGravitational3D1(types_coords,types_dict)
            res['rygr_'+atom_type]=CalculateRadiusofGyration(types_coords,types_dict)
            res['Harary3D_'+atom_type]=CalculateHarary3D(DistanceMatrix)
            res['AGDD_'+atom_type]=CalculateAverageGeometricalDistanceDegree(DistanceMatrix)
            res['SEig_'+atom_type]=CalculateAbsEigenvalueSumOnGeometricMatrix(DistanceMatrix)
            res['SPAN_'+atom_type]=CalculateSPANR(types_coords,types_dict)
            res['ASPAN_'+atom_type]=CalculateAverageSPANR(types_coords,types_dict)
            res['MEcc_'+atom_type]=CalculateMolecularEccentricity(types_coords,types_dict)
#            res['IA_'+atom_type]=CalculatePrincipalMomentofInertia(types_coords,types_dict)['IA']
#            res['IB_'+atom_type]=CalculatePrincipalMomentofInertia(types_coords,types_dict)['IB']
#            res['IC_'+atom_type]=CalculatePrincipalMomentofInertia(types_coords,types_dict)['IC']
            res['IA/B_'+atom_type]=CalculateRatioPMI(types_coords,types_dict)['IA/B']
            res['IA/C_'+atom_type]=CalculateRatioPMI(types_coords,types_dict)['IA/C']
            res['IB/C_'+atom_type]=CalculateRatioPMI(types_coords,types_dict)['IB/C']
        elif not types_coords:
            res={}
            res['W3D_'+atom_type]=0
            res['Petitj3D_'+atom_type]=0
            res['GeDi_'+atom_type]=0
            res['grav_'+atom_type]=0
            res['rygr_'+atom_type]=0
            res['Harary3D_'+atom_type]=0
            res['AGDD_'+atom_type]=0
            res['SEig_'+atom_type]=0
            res['SPAN_'+atom_type]=0
            res['ASPAN_'+atom_type]=0
            res['MEcc_'+atom_type]=0
#            res['IA_'+atom_type]=0
#            res['IB_'+atom_type]=0
#            res['IC_'+atom_type]=0
            res['IA/B_'+atom_type]=0
            res['IA/C_'+atom_type]=0
            res['IB/C_'+atom_type]=0
    elif atom_type in ('no_type'):
        DistanceMatrix=GetGementricalDistanceMatrix(Coordinates)
        res={}
        res['W3D_'+atom_type]=Calculate3DWiener(DistanceMatrix)
        res['Petitj3D_'+atom_type]=CalculatePetitjean3DIndex(DistanceMatrix)
        res['GeDi_'+atom_type]=CalculateGemetricalDiameter(DistanceMatrix)
        res['grav_'+atom_type]=CalculateGravitational3D1(Coordinates,types_dict)
        res['rygr_'+atom_type]=CalculateRadiusofGyration(Coordinates,types_dict)
        res['Harary3D_'+atom_type]=CalculateHarary3D(DistanceMatrix)
        res['AGDD_'+atom_type]=CalculateAverageGeometricalDistanceDegree(DistanceMatrix)
        res['SEig_'+atom_type]=CalculateAbsEigenvalueSumOnGeometricMatrix(DistanceMatrix)
        res['SPAN_'+atom_type]=CalculateSPANR(Coordinates,types_dict)
        res['ASPAN_'+atom_type]=CalculateAverageSPANR(Coordinates,types_dict)
        res['MEcc_'+atom_type]=CalculateMolecularEccentricity(Coordinates,types_dict)
#        res['IA_'+atom_type]=CalculatePrincipalMomentofInertia(Coordinates,types_dict)['IA']
#        res['IB_'+atom_type]=CalculatePrincipalMomentofInertia(Coordinates,types_dict)['IB']
#        res['IC_'+atom_type]=CalculatePrincipalMomentofInertia(Coordinates,types_dict)['IC']
        res['IA/B_'+atom_type]=CalculateRatioPMI(Coordinates,types_dict)['IA/B']
        res['IA/C_'+atom_type]=CalculateRatioPMI(Coordinates,types_dict)['IA/C']
        res['IB/C_'+atom_type]=CalculateRatioPMI(Coordinates,types_dict)['IB/C']
    return res

def GetGeometric_alpha_coarse(Coordinates,types_dict):
    result={}
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'no_type'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'HP'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'PL'))
    return result

def GetGeometric_alpha_lig_coarse(Coordinates,types_dict):
    result={}
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'no_type'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'HP'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'PL'))
    return result

def GetGeometric_alpha_fine(Coordinates,types_dict):
    result={}
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'no_type'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'AR'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'HP'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'DA'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'A'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'D'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'PL'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'P'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'N'))
    return result

def GetGeometric_alpha_lig_fine(Coordinates,types_dict):
    result={}
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'no_type'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'AR'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'HP'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'DA'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'A'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'D'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'PL'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'P'))
    result.update(GetGeometric_withtypes(Coordinates,types_dict,'N'))
#    result.update(GetGeometric_withtypes(Coordinates,types_dict,'HX'))
    return result

if __name__=="__main__":

    print(res)


