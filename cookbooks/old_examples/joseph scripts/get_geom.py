# -*- coding: utf-8 -*-

import os 
import glob
import itertools
import numpy as np
import numpy
import scipy
import math
import scipy.cluster.hierarchy as hier
from operator import itemgetter
from collections import defaultdict
from collections import OrderedDict
from numpy import linalg, mean, tile, mat, transpose

## This gets the alpha atom coordinates##
def get_alpha_coordinates(pdb_directory):
    Coordinates=[]
    if pdb_directory.endswith(".pdb"):
        with open(pdb_directory,'r') as f:
            templines=f.readlines()
        for i in templines:
            if i[17:20]=='AAC':
                Coordinates.append([i[17:20].strip(),i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
            else:
                pass
    else:
        for filename in sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb')):
            with open(filename,'r') as f:
                templines=f.readlines()
            for i in templines:
                if i[17:20]=='AAC':
                    Coordinates.append([i[17:20].strip(),i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
                else:
                    pass
        Coordinates.sort()
        Coordinates=list(Coordinates for Coordinates,_ in itertools.groupby(Coordinates))
    return Coordinates



#def protein_parser(filename):
#    parsed_data = defaultdict(dict)
#    with open(filename, 'r') as f:
#        templines= f.readlines()
#    for i in templines:
#        res_id=i[22:26].strip()
#        atom_name=i[13:16].strip()
#        parsed_data[res_id][atom_name] = [i[0:6].strip(),i[6:11].strip(),i[13:16].strip(),i[17:20].strip(),i[21].strip(),i[22:26].strip(),i[30:38].strip(),i[38:46].strip(),i[46:54].strip(),i[54:62].strip()]
#    return parsed_data

#    with open(filename, 'r') as f:
#        lines = f.readlines()
#    for line in lines:
#        line = line.split()
#        parsed_data[line[5]][line[2]] =[
#    return parsed_data

def get_surf_from_pdb(pdb_directory):

    Coordinates=[]
    if pdb_directory.endswith(".pdb"):
        with open(pdb_directory,'r') as f:
            templines=f.readlines()
        for i in templines:
            temp=i.split()
            if i[17:20]=='AAC':
                break
            else:
                Coordinates.append([i[17:20].strip(),i[12:16].strip(),i[61:].strip()])
    else:
        file_list=sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb'))
        for i,file in enumerate(file_list):
            if file.endswith("beta.pdb"):
                file_list.pop(i)
        file_list.pop(-1)
        for filename in file_list:
            with open(filename,'r') as f:
                templines=f.readlines()
            for i in templines:
    #            temp=i.split()
                if i[17:20]=='AAC':
                    break
                else:
                    Coordinates.append([i[17:20].strip(),i[12:16].strip(),i[61:].strip(),int(i[7:12].strip())])
    #                Coordinates.append([temp[3],temp[5],temp[2]])    
        Coordinates.sort()
        toret=OrderedDict()
        for i in Coordinates:
            if i[3] not in toret:
                toret[i[3]]=i
        Coordinates=[toret[keys][0:3] for keys in toret.keys()]
    return Coordinates

## this gets scores from the pdb 

def get_score_from_pdb(pdb_directory):

    Coordinates=[]
    if pdb_directory.endswith(".pdb"):
        with open(pdb_directory,'r') as f:
            templines=f.readlines()
        for i in templines:
            temp=i.split()
            if i[17:20]=='AAC':
                break
            else:
                Coordinates.append([i[17:20].strip(),i[12:16].strip(),i[55:60].strip()])
    else:
        file_list=sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb'))
        for i,file in enumerate(file_list):
            if file.endswith("beta.pdb"):
                file_list.pop(i)
        file_list.pop(-1)
        for filename in file_list:
            with open(filename,'r') as f:
                templines=f.readlines()
            for i in templines:
    #            temp=i.split()
                if i[17:20]=='AAC':
                    break
                else:
                    Coordinates.append([i[17:20].strip(),i[12:16].strip(),i[55:60].strip(),int(i[7:12].strip())])
    #                Coordinates.append([temp[3],temp[5],temp[2]])    
        Coordinates.sort()
        toret=OrderedDict()
        for i in Coordinates:
            if i[3] not in toret:
                toret[i[3]]=i
        Coordinates=[toret[keys][0:3] for keys in toret.keys()]
    return Coordinates



def get_coordinates_from_pdb(pdb_directory):
    Coordinates=[]
    if pdb_directory.endswith(".pdb"):
        with open(pdb_directory,'r') as f:
            templines=f.readlines()
        for i in templines:
            temp=i.split()
            if i[17:20]=='AAC':
                break
            else:
                Coordinates.append([i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
    else:
        file_list=sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb'))
        for i,file in enumerate(file_list):
            if file.endswith("beta.pdb"):
                file_list.pop(i)
        file_list.pop(-1)
        for filename in file_list:
            with open(filename,'r') as f:
                templines=f.readlines()
            for i in templines:
                temp=i.split()
                if i[17:20]=='AAC':
                    break
                else:
                    Coordinates.append([i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
        Coordinates.sort()
        Coordinates=list(Coordinates for Coordinates,_ in itertools.groupby(Coordinates))
    return Coordinates

## This gets gets the Coordinates file with the atom type from the ac file
def get_lig_coordinates_from_ac(pdb_file):
    os.system("obabel " + pdb_file + " -O temp.pdb -p 7")
    os.system("~/amber14/bin/antechamber -fi pdb -i temp.pdb -o lig.ac -fo ac")
    with open('lig.ac','r') as f:
        templines=f.readlines()
    Coordinates=[]
    templines=templines[2:]
    for i in templines:
        if i[12:16].strip()=='H':
            break
        elif i[0:4]=='BOND':
            break
        else:
            Coordinates.append(['LIG',i[72:].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
    os.system('rm ANTECHAMBER* ATOMTYPE* lig.ac temp.pdb')
    return Coordinates

## get ligands surface from nacces calc
def get_lig_surf_from_ac(pdb_file):
    os.system("/Users/jkatigba/bin/Naccess/naccess -h " + pdb_file)
    os.system("obabel " + pdb_file + " -O temp.pdb -p 7")
    os.system("~/amber14/bin/antechamber -fi pdb -i temp.pdb -o lig.ac -fo ac")
    with open('lig.ac','r') as f:
        ac_temp=f.readlines()
    with open('lig.asa','r') as f:
        asa_temp=f.readlines()
    Coordinates=[]
    ac_temp=ac_temp[2:]
    for index,i in enumerate(ac_temp):
        if i[12:16].strip()=='H':
            break
        elif i[0:4]=='BOND':
            break
        else:
            Coordinates.append(['LIG',i[72:].strip(),asa_temp[index][56:63].strip()])
    os.system('rm ANTECHAMBER* ATOMTYPE* lig.ac temp.pdb lig.log lig.asa lig.rsa')
    return Coordinates

def GetAtomDistance(x,y):
    temp=[math.pow(x[0]-y[0],2),math.pow(x[1]-y[1],2),math.pow(x[2]-y[2],2)]
    res=scipy.sqrt(sum(temp))
    return res


def _GetGeometricalCenter(ChargeCoordinates):
    res1=[]
    res2=[]
    res3=[]
    for i in ChargeCoordinates:
        res1.append(float(i[2]))
        res2.append(float(i[3]))
        res3.append(float(i[4]))
    result=[scipy.mean(res1),scipy.mean(res2),scipy.mean(res3)]
    return result


def get_types_dict(filename):
    types_dict = defaultdict(dict)
    fp = open(filename, 'r')
    lines = fp.readlines()
    fp.close()
    for line in lines:
        line = line.split()
        types_dict[line[2]][line[0]]=line

    return types_dict

def GetGementricalDistanceMatrix(Coordinates):

    Coords=[]
    for i in Coordinates:
        Coords.append([float(i[2]),float(i[3]),float(i[4])])
    NAtom=len(Coords)
    DistanceMatrix=scipy.zeros((NAtom,NAtom))
    for i in range(NAtom-1):
        for j in range(i+1,NAtom):
            DistanceMatrix[i,j]=GetAtomDistance(Coords[i],Coords[j])
            DistanceMatrix[j,i]=DistanceMatrix[i,j]
    return DistanceMatrix

def get_all_res_coordinates(Coordinates,parsed_data):
    all_res=[]
    for i in Coordinates:
        for j in parsed_data[i[1]].keys():
#            all_res.append([parsed_data[i[1]][j][17:20],parsed_data[i[1]][j][22:26].strip(),parsed_data[i[1]][j][12:16].strip(),parsed_data[i[1]][j][31:38].strip(),parsed_data[i[1]][j][39:46].strip(),parsed_data[i[1]][j][47:54].strip()])
            all_res.append([parsed_data[i[1]][j][3],parsed_data[i[1]][j][2],parsed_data[i[1]][j][6],parsed_data[i[1]][j][7],parsed_data[i[1]][j][8]])

    all_res.sort()
    all_res=list(all_res for all_res,_ in itertools.groupby(all_res))

    return all_res

def get_coordinates_from_pdb_resid(pdb_directory):
    Coordinates=[]
    for filename in sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb')):
        with open(filename,'r') as f:
            templines=f.readlines()
        for i in templines:
            if i[17:20]=='AAC':
                break
            else:
                Coordinates.append([i[17:20],i[22:26].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
    Coordinates.sort()
    Coordinates=list(Coordinates for Coordinates,_ in itertools.groupby(Coordinates))
    return Coordinates

def get_coordinates_bypock_resid(filename):
    Coordinates=[]
    with open(filename,'r') as f:
        templines=f.readlines()
    for i in templines:
        if i[17:20]=='AAC':
            break
        else:
            Coordinates.append([i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),i[22:26].strip()])
    return Coordinates




## get_type_coords
# coarse grained types Hydrophobic, Hydrophilic
# fine grained types Charge+,Charge-,hbd,hba,doneptor,polar,aromatic,aliphatic
def get_type_coords(Coordinates,types_dict,atom_types):
    types_coords=[]
    for i in Coordinates:
        if types_dict[i[0]][i[1]][4]==atom_types:
            types_coords.append(i)    
    return types_coords

## get type surfaces

def get_type_surf(Surface,types_dict,atom_types):
    types_surf=[]
    for i in Surface:
        if types_dict[i[0]][i[1]][4]==atom_types:
            types_surf.append(i)
    return types_surf

## this gets type_scores
def get_type_score(Score,types_dict,atom_types):
    types_score=[]
    for i in Score:
        if types_dict[i[0]][i[1]][4]==atom_types:
            types_score.append(i)
    return types_score

## change the thing
# filename is table_pockets.dat
def geom_screen_by_occupancy(filename,directory,occupancy):
    with open(filename,'r') as f:
        templines=f.readlines()
    pock_list=[]
    templines.pop(0)
    for i in templines:
        temp=i.split()
        if float(temp[3].strip('%'))/100 > float(occupancy):
            pock_list.append('pockets/'+str(int(temp[0])-1).zfill(3)+'.pdb')
    Coords=[]
    for j in pock_list:
        Coords.extend(get_coordinates_from_pdb(directory.strip('/')+'/'+j))
    Coords.sort()
    Coords=list(Coords for Coords,_ in itertools.groupby(Coords))
    return Coords

def surf_screen_by_occupancy(filename,directory,occupancy):
    with open(filename,'r') as f:
        templines=f.readlines()
    pock_list=[]
    templines.pop(0)
    for i in templines:
        temp=i.split()
        if float(temp[3].strip('%'))/100 > float(occupancy):
            pock_list.append(directory+'/pockets/'+str(int(temp[0])-1).zfill(3)+'.pdb')
    Coords=[]

    for pockname in pock_list:
        with open(pockname,'r') as f:
            templines=f.readlines()
        for i in templines:
            if i[17:20]=='AAC':
                break
            else:
                Coords.append([i[17:20].strip(),i[12:16].strip(),i[61:].strip(),int(i[7:12].strip())])
    Coords.sort()
    toret=OrderedDict()
    for i in Coords:
        if i[3] not in toret:
            toret[i[3]]=i
    Coords=[toret[keys][0:3] for keys in toret.keys()]
    return Coords

def score_screen_by_occupancy(filename,directory,occupancy):
    with open(filename,'r') as f:
        templines=f.readlines()
    pock_list=[]
    templines.pop(0)
    for i in templines:
        temp=i.split()
        if float(temp[3].strip('%'))/100 > float(occupancy):
            pock_list.append(directory+'/pockets/'+str(int(temp[0])-1).zfill(3)+'.pdb')
    Coords=[]

    for pockname in pock_list:
        with open(pockname,'r') as f:
            templines=f.readlines()
        for i in templines:
            if i[17:20]=='AAC':
                break
            else:
                Coords.append([i[17:20].strip(),i[12:16].strip(),i[55:60].strip(),int(i[7:12].strip())])
    Coords.sort()
    toret=OrderedDict()
    for i in Coords:
        if i[3] not in toret:
            toret[i[3]]=i
    Coords=[toret[keys][0:3] for keys in toret.keys()]
    return Coords


## ligand_fragmentizer
## submit list of fragments names
#lig_list=[['N6','C6','C5','N7','C8','N9','C4','N3','C2','N1'],["C1'","O4'","C4'","C5'","C3'","O3'","C2'","O2'"],["O5'","PA","PB","PG","O1A","O2A","O3A","O1B","O2B","O3B","O1G","O2G","O3G"]]
#lig_name='ATP'
#get_LIG_subdivide(filename,directory,lig_list,lig_name)


def get_LIG_subdivide(filename,directory,lig_list,lig_name):
    with open(filename,'r') as f:
        templines=f.readlines()
    prot=[]
    for i in templines:
        if i[0:6]=='ATOM  ':
            prot.append(i)
    prot.append('TER \n')
    lig_dict={}
    cwd=directory.rstrip('/')
    
    for ix,i in enumerate(lig_list):
        temp=[]
        for j in i:
            temp.append(next(k for k in templines if k[12:16].strip()==j and k[17:20]==lig_name))
        lig_dict['fragment_'+str(ix)]=temp
        os.mkdir(cwd+'/fragment_'+str(ix))
        temp_lig_prot=[]
        temp_lig_prot.extend(prot)
        temp_lig_prot.extend(temp)
        with open(cwd+'/fragment_'+str(ix)+'/fragment_'+str(ix)+'.pdb','w') as g:
            g.writelines(temp_lig_prot)

# gets coordinates for specific pocket and fuses them together
def fuse_coord_list(pock_list):
    coords=[]
    for i in pock_list:
        coords.extend(get_coordinates_from_pdb(i+'.pdb'))
    coords.sort()
    coords=list(coords for coords,_ in itertools.groupby(coords))
    return coords


# soergel distance
def soergel(i,j):
    li=len(i)
    lj=len(j)
    if li is not lj:
        print('lengths not equal')
    else:
        ntemp=0
        dtemp=0
        for k in range(li):
            ntemp=ntemp+abs(i[k]-j[k])
            dtemp=dtemp+max(i[k],j[k])
        score=ntemp/dtemp
    return score
    
def fuse_lig_beta(lig_file,beta_file,file_name):
    with open(lig_file,'r') as f:
        templines=f.readlines()
    templines.append("TER\n")
    with open(beta_file,'r') as f:
        betalines=f.readlines()
    templines.extend(betalines)
    with open(file_name,'w') as g:
        g.writelines(templines)

def write_beta_only(beta_file,file_name):
    with open(beta_file,'r') as f:
        betalines=f.readlines()
    Coordinates=[]
    for i in betalines:
        if i[13:16]=='BAO':
            Coordinates.append(i[:54]+'                       C \n')
    with open(file_name,'w') as g:
        g.writelines(Coordinates)

def write_alpha_only(pdb_directory,file_name):
    Coordinates=[]
    file_list=sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb'))
    for i,file in enumerate(file_list):
        if file.endswith("beta.pdb"):
            file_list.pop(i)
    file_list.pop(-1)
    for filename in file_list:
        with open(filename,'r') as f:
            templines=f.readlines()
        for i in templines:
            temp=i.split()
            if i[17:20]=='AAC':
                Coordinates.append(i[:54]+'                       C \n')
    with open(file_name,'w') as g:
        g.writelines(Coordinates)

def write_beta_ACC(beta_file,pdb_directory,file_name):
    with open(beta_file,'r') as f:
        betalines=f.readlines()
    Coordinates=[]
    for i in betalines:
        if i[17:20]=='BAC':
            Coordinates.append(i[:54]+'                       C \n')
    file_list=sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb'))
    for i,file in enumerate(file_list):
        if file.endswith("beta.pdb"):
            file_list.pop(i)
    file_list.pop(-1)
    for filename in file_list:
        with open(filename,'r') as f:
            templines=f.readlines()
        for i in templines:
            if i[17:20]=='ACC':
                Coordinates.append(i[:54]+'                       B \n')
    with open(file_name,'w') as g:
        g.writelines(Coordinates)

def fuse_lig_beta_bcc(lig_file,beta_file_list,filename):
    with open(lig_file,'r') as f:
        templines=f.readlines()
    templines.append("TER\n")
    
    beta_templines=[]
    for beta_name in beta_file_list:
        with open(beta_name,'r') as f:
            betalines=f.readlines()
        for i in betalines:
            if i[17:20]=='BCC':
                beta_templines.append(i)
    
    beta_temp=[]
    for ix,i in enumerate(beta_templines):
        beta_temp.append(i[0:22]+str(ix).rjust(4)+i[26:])
    
    templines.extend(beta_temp)
    
    with open(file_name,'w') as g:
        g.writelines(templines)


def transform_lig(trans_coords,untrans_coords,lig_file,out_file):
    A=[]
    for i in untrans_coords:
        A.append([float(i[2]),float(i[3]),float(i[4])])
    B=[]
    for j in trans_coords:
        B.append([float(j[2]),float(j[3]),float(j[4])])
    A=mat(A) 
    B=mat(B)
    assert len(A) == len(B)

    N = A.shape[0]; # total points

    centroid_A = mean(A, axis=0)
    centroid_B = mean(B, axis=0)

    # centre the points
    AA = A - tile(centroid_A, (N, 1))
    BB = B - tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = transpose(AA) * BB

    U, S, Vt = linalg.svd(H)

    R = Vt.T * U.T

    # special reflection case
    if linalg.det(R) < 0:
       print("Reflection detected")
       Vt[2,:] *= -1
       R = Vt.T * U.T

    t = -R*centroid_A.T + centroid_B.T
    lig_coords=get_lig_coordinates_from_ac(lig_file)
    C=[]
    for k in lig_coords:
        C.append([float(k[2]),float(k[3]),float(k[4])])
    D=(R*mat(C).T)+tile(t,(1,len(C)))
    D=D.T
    x=len(D)
    with open(lig_file,'r') as f:
        templines=f.readlines() 
    transform=[]
    for lx,l in enumerate(templines):
        transform.append(l[:30]+str(round(D[lx,0],3)).rjust(8)+str(round(D[lx,1],3)).rjust(8)+str(round(D[lx,2],3)).rjust(8)+l[54:])
    with open(out_file,'w') as g:
        g.writelines(transform) 

def get_aligned_beta_coords(pdb_file):
    with open(pdb_file,"r") as f:
        templines=f.readlines()
    Coordinates=[]
    templines.pop(0)
    templines.pop(0)
    for i in templines:
        if i[0:6]=='HETATM':
            Coordinates.append([i[17:20],i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
    return Coordinates

def get_beta_coord_list_all(beta_pdb):
    with open(beta_pdb,"r") as f:
        templines=f.readlines()
    Coordinates=[]
    ix=0.0
    for i in templines:
        if i[17:20]=='BAC':
#            if i[25]=='1':
            Coordinates.append([i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),float(i[66:72].strip()),str(int(ix))])
            ix=ix+1
#            Coordinates.append([i[17:20],i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),"Hydrophobic",i[55:60].strip(),i[61:64].strip()])
#            elif i[25]=='2':
#                Coordinates.append([i[17:20],i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),"Polar",i[55:60].strip(),i[61:64].strip()])
#            elif i[25]=='3':
#                Coordinates.append([i[17:20],i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),"Charge",i[55:60].strip(),i[61:64].strip()])

    return Coordinates

def get_beta_coord_list_occ(beta_pdb):
    with open(beta_pdb,"r") as f:
        templines=f.readlines()
    Coordinates=[]
    ix=0.0
    for i in templines:
        if i[13:16]=='BAO':
            Coordinates.append([i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),float(i[66:72].strip()),str(int(ix))])
            ix=ix+1
#            if i[25]=='1':
#                Coordinates.append([i[17:20],i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),"Hydrophobic",i[55:60].strip(),i[61:64].strip()])
#            elif i[25]=='2':
#                Coordinates.append([i[17:20],i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),"Polar",i[55:60].strip(),i[61:64].strip()])
#            elif i[25]=='3':
#                Coordinates.append([i[17:20],i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),"Charge",i[55:60].strip(),i[61:64].strip()])

    return Coordinates

def get_ACC_PAC(pdb_directory):
    PAC_Coords=[]
    ACC_Coords=[]
    Coordinates=[]
    if pdb_directory.endswith(".pdb"):
        with open(pdb_directory,'r') as f:
            templines=f.readlines()
        for i in templines:
            if i[17:20]=='ACC':
                ACC_Coords.append([float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
            elif i[17:20]=='PAC':
                PAC_Coords.append([float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
        Coordinates=[ACC_Coords,PAC_Coords]
    else:
        file_list=sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb'))
        for i,file in enumerate(file_list):
            if file.endswith("beta.pdb"):
                file_list.pop(i)
        file_list.pop(-1)
        for filename in file_list:
            with open(filename,'r') as f:
                templines=f.readlines()
            for i in templines:
                if i[17:20]=='ACC':
                    ACC_Coords.append([float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
                elif i[17:20]=='PAC':
                    PAC_Coords.append([float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
        for i in range(len(file_list)):
            Coordinates.append([ACC_Coords[i],PAC_Coords[i]])
#    Coordinates=[x for x in Coordinates if x != []]
    return Coordinates


# get coordinates with res id and atom name

def get_coordinates_from_pdb_withID(pdb_directory):
    Coordinates=[]
    if pdb_directory.endswith(".pdb"):
        with open(pdb_directory,'r') as f:
            templines=f.readlines()
        for i in templines:
            temp=i.split()
            if i[17:20]=='AAC':
                break
            else:
                Coordinates.append([i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),i[22:26].strip()])
    else:
        file_list=sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb'))
        for i,file in enumerate(file_list):
            if file.endswith("beta.pdb"):
                file_list.pop(i)
        file_list.pop(-1)
        for filename in file_list:
            with open(filename,'r') as f:
                templines=f.readlines()
            for i in templines:
                temp=i.split()
                if i[17:20]=='AAC':
                    break
                else:
                    Coordinates.append([i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),i[22:26].strip()])
        Coordinates.sort()
        Coordinates=list(Coordinates for Coordinates,_ in itertools.groupby(Coordinates))
    return Coordinates

# get parsed_dict

def protein_sasa_parser(filename):
    parsed_data = defaultdict(dict)
    with open(filename, 'r') as f:
        templines= f.readlines()
    for i in templines:
        if i[0:6]=='ATOM  ' or i[0:6]=='HETATM':
            res_id=i[22:26].strip()
            res_name=i[17:20].strip()
            atom_name=i[12:16].strip()
            parsed_data[res_id][atom_name] = [i[17:20].strip(),i[12:16].strip(),float(i[54:62].strip()),i[22:26].strip()]
        elif i[0:6]=='CONECT':
            break
    return parsed_data

def protein_parser(filename):
    parsed_data = defaultdict(dict)
    with open(filename, 'r') as f:
        templines= f.readlines()
    for i in templines:
        if i[0:6]=='ATOM  ' or i[0:6]=='HETATM':
            res_id=i[22:26].strip()
            res_name=i[17:20].strip()
            atom_name=i[12:16].strip()
            parsed_data[res_id][atom_name] = [i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),i[22:26].strip()]
        elif i[0:6]=='CONECT':
            break
    return parsed_data

# this makes updates the protonation and resname in dict
def update_prot_coords_name(prot_coords,parsed_dict):
    for ix,i in enumerate(prot_coords):
        if i[0] in ['CYS','HIS']:
            new_name=parsed_dict[i[-1]]['C'][0]
            prot_coords[ix][0]=new_name


# this updates the coordslist to include carbonyl backbone carbons

def add_c_coords(prot_coords,parsed_dict):
    mappings={'ALA':{'O':['C']},
              'GLY':{'O':['C']},
              'SER':{'O':['C']},
              'THR':{'O':['C']},
              'LEU':{'O':['C']},
              'ILE':{'O':['C']},
              'VAL':{'O':['C']},
              'ASN':{'O':['C']},
              'GLN':{'O':['C']},
              'ARG':{'O':['C']},
              'HID':{'O':['C']},
              'HIE':{'O':['C']},
              'HIS':{'O':['C']},
              'HIP':{'O':['C']},
              'TRP':{'O':['C']},
              'PHE':{'O':['C']},
              'TYR':{'O':['C']},
              'GLU':{'O':['C']},
              'ASP':{'O':['C']},
              'LYS':{'O':['C']},
              'LYN':{'O':['C']},
              'CYS':{'O':['C']},
              'CYM':{'O':['C']},
              'CYX':{'O':['C']},
              'MET':{'O':['C']},
              'PRO':{'O':['C']},
              'ASH':{'O':['C']},
              'GLH':{'O':['C']}
              }


    c_coords=[]
    for i in prot_coords:
        try:
            for j in mappings[i[0]][i[1]]:
                c_coords.append(parsed_dict[i[-1]][j])
        except KeyError:
            continue

    prot_coords.extend(c_coords)

    prot_coords.sort()
    prot_coords=list(prot_coords for prot_coords,_ in itertools.groupby(prot_coords))
    return prot_coords


# this updates the coords list to include hydrogens
def add_h_coords(prot_coords,parsed_dict):
    mappings={'ALA':{'N':['H']},
              'GLY':{'N':['H']},
              'SER':{'N':['H'],'OG':['HG']},
              'THR':{'N':['H'],'OG1':['HG1']},
              'LEU':{'N':['H'],},
              'ILE':{'N':['H'],},
              'VAL':{'N':['H'],},
              'ASN':{'N':['H'],'ND2':['HD21','HD22']},
              'GLN':{'N':['H'],'NE2':['HE21','HE22']},
              'ARG':{'N':['H'],'NE':['HE'],'NH1':['HH11','HH12'],'NH2':['HH21','HH22']},
              'HID':{'N':['H'],'ND1':['HD1']},
              'HIE':{'N':['H'],'NE2':['HE2']},
              'HIS':{'N':['H'],'ND1':['HD1'],'NE2':['HE2']},
              'HIP':{'N':['H'],'ND1':['HD1'],'NE2':['HE2']},
              'TRP':{'N':['H'],'NE1':['HE1']},
              'PHE':{'N':['H']},
              'TYR':{'N':['H'],'OH':['HH']},
              'GLU':{'N':['H']},
              'ASP':{'N':['H']},
              'LYS':{'N':['H'],'NZ':['HZ1','HZ2','HZ3']},
              'LYN':{'N':['H'],'NZ':['HZ2','HZ3']},
              'CYS':{'N':['H'],'SG':['HG']},
              'CYM':{'N':['H']},
              'CYX':{'N':['H']},
              'MET':{'N':['H']},
              'ASH':{'N':['H'],'OD2':['HD2']},
              'GLH':{'N':['H'],'OE2':['HE2']}
              }

#    mappings={'ALA':{'N':['H'],'CA':['HA'],'CB':['HB1','HB2','HB3']},
#              'GLY':{'N':['H'],'CA':['HA2','HA3']},
#              'SER':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'OG':['HG']},
#              'THR':{'N':['H'],'CA':['HA'],'CB':['HB'],'OG1':['HG1'],'CG2':['HG21','HG22','HG23']},
#              'LEU':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CG':['HG'],'CD1':['HD11','HD12','HD13'],'CD2':['HD21','HD22','HD23']},
#              'ILE':{'N':['H'],'CA':['HA'],'CB':['HB'],'CG2':['HG21','HG22','HG23'],'CG1':['HG12','HG13'],'CD1':['HD11','HD12','HD13']},
#              'VAL':{'N':['H'],'CA':['HA'],'CB':['HB'],'CG2':['HG21','HG22','HG23'],'CG1':['HG11','HG12','HG13']},
#              'ASN':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'ND2':['HD21','HD22']},
#              'GLN':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CG':['HG2','HG3'],'NE2':['HE21','HE22']},
#              'ARG':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CG':['HG2','HG3'],'CD':['HD2','HD3'],'NE':['HE'],'NH1':['HH11','HH12'],'NH2':['HH21','HH22']},
#              'HID':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'ND1':['HD1'],'CE1':['HE1'],'CD2':['HD2']},
#              'HIE':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'NE2':['HE2'],'CE1':['HE1'],'CD2':['HD2']},
#              'HIS':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'ND1':['HD1'],'CE1':['HE1'],'CD2':['HD2']},
#              'HIP':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'ND1':['HD1'],'NE2':['HE2'],'CE1':['HE1'],'CD2':['HD2']},
#              'TRP':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CD1':['HD1'],'NE1':['HE1'],'CZ2':['HZ2'],'CH2':['HH2'],'CZ3':['HZ3'],'CE3':['HE3']},
#              'PHE':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CD1':['HD1'],'CE1':['HE1'],'CZ':['HZ'],'CE2':['HE2'],'CD2':['HD2']},
#              'TYR':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CD1':['HD1'],'CE1':['HE1'],'OH':['HH'],'CE2':['HE2'],'CD2':['HD2']},
#              'GLU':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CG':['HG2','HG3']},
#              'ASP':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3']},
#              'LYS':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CG':['HG2','HG3'],'CD':['HD2','HD3'],'CE':['HE2','HE3'],'NZ':['HZ1','HZ2','HZ3']},
#              'LYN':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CG':['HG2','HG3'],'CD':['HD2','HD3'],'CE':['HE2','HE3'],'NZ':['HZ2','HZ3']},
#              'PRO':{'CA':['HA'],'CB':['HB2','HB3'],'CG':['HG2','HG3'],'CD':['HD2','HD3']},
#              'CYS':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'SG':['HG']},
#              'CYM':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3']},
#              'CYX':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3']},
#              'MET':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CG':['HG2','HG3'],'CE':['HE1','HE2','HE3']},
#              'ASH':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'OD2':['HD2']},
#              'GLH':{'N':['H'],'CA':['HA'],'CB':['HB2','HB3'],'CG':['HG2','HG3'],'OE2':['HE2']}
#              }

    h_coords=[]
    for i in prot_coords:
        try: 
            for j in mappings[i[0]][i[1]]:
                h_coords.append(parsed_dict[i[-1]][j])
        except KeyError:
            continue
    
    prot_coords.extend(h_coords)        
    return prot_coords

        
def checking_free_polar_sites_withC(prot_coords,parsed_dict):
    hba_list=['O']
    hbd_list=['N']
    res_list=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','HID','HIE','HIP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','ASH','GLH','CYM','CYX','LYN']
    

    for ix,i in enumerate(prot_coords):
        if i[1] in hba_list and i[0] in res_list:
            i_coords=np.array([float(i[2]),float(i[3]),float(i[4])])
            for j in parsed_dict.keys():
                if 'N' in parsed_dict[j].keys() and parsed_dict[j]['N'][0] in res_list:
                    j_coords=np.array([float(parsed_dict[j]['N'][2]),float(parsed_dict[j]['N'][3]),float(parsed_dict[j]['N'][4])])
                    if GetAtomDistance(j_coords,i_coords) < 3.4 and parsed_dict[j]['N'][-1] != i[-1]:
                        c_coords=np.array([float(parsed_dict[i[-1]]['C'][2]),float(parsed_dict[i[-1]]['C'][3]),float(parsed_dict[i[-1]]['C'][4])])
                        ci=c_coords-i_coords
                        ji=j_coords-i_coords
                        cosine_angle = np.dot(ci, ji) / (np.linalg.norm(ci) * np.linalg.norm(ji))
                        angle = np.arccos(cosine_angle)*57.2958
                        if angle >= 90:
                            prot_coords[ix][1]='OP'
    #                    print i,parsed_dict[j]['N'],GetAtomDistance(j_coords,i_coords),angle
        elif i[1] in hbd_list and i[0] in res_list:
            i_coords=np.array([float(i[2]),float(i[3]),float(i[4])])
            for j in parsed_dict.keys():
                if 'O' in parsed_dict[j].keys() and parsed_dict[j]['O'][0] in res_list:    
                    j_coords=np.array([float(parsed_dict[j]['O'][2]),float(parsed_dict[j]['O'][3]),float(parsed_dict[j]['O'][4])])
                    if GetAtomDistance(j_coords,i_coords) < 3.4 and parsed_dict[j]['O'][-1] != i[-1]:
                        c_coords=np.array([float(parsed_dict[j]['C'][2]),float(parsed_dict[j]['C'][3]),float(parsed_dict[j]['C'][4])])
                        ci=c_coords-j_coords
                        ji=i_coords-j_coords
                        cosine_angle = np.dot(ci, ji) / (np.linalg.norm(ci) * np.linalg.norm(ji))
                        angle = np.arccos(cosine_angle)*57.2958
                        if angle >= 90:
                            prot_coords[ix][1]='NP'
    #                    print i,parsed_dict[j]['O'],GetAtomDistance(j_coords,i_coords),angle

    return prot_coords

#gets angle, must feed p1,p2,c in np.array. p1,p2 are points and c is center
def get_angle(p1,p2,c):
    p1c=p1-c
    p2c=p2-c
    cosine_angle = np.dot(p1c, p2c) / (np.linalg.norm(p1c) * np.linalg.norm(p2c))
    angle = np.arccos(cosine_angle)*57.2958
    return angle

# gets dihedral angle, p1 and p2 are the central axis, p0 and p3 are the end points
def get_dihedral(p0,p1,p2,p3):

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

## check if beta atom has the proper angle for an H-bond, only works for amino acids. Must filter first to make sure only coords that are amino acids

def beta_prot_angle_check(p_coord,b_coord,prot_dict,parsed_dict):
  
    A_mappings={'ALA':{'O':['C'],'OXT':['C']},
              'GLY':{'O':['C'],'OXT':['C']},
              'SER':{'O':['C'],'OXT':['C'],'OG':['CB']},
              'THR':{'O':['C'],'OXT':['C'],'OG1':['CB']},
              'LEU':{'O':['C'],'OXT':['C']},
              'ILE':{'O':['C'],'OXT':['C']},
              'VAL':{'O':['C'],'OXT':['C']},
              'ASN':{'O':['C'],'OXT':['C'],'OD1':['CG']},
              'GLN':{'O':['C'],'OXT':['C'],'OE1':['CD']},
              'ARG':{'O':['C'],'OXT':['C']},
              'HID':{'O':['C'],'OXT':['C'],'NE2':['CE1'],'ND1':['CG']},
              'HIE':{'O':['C'],'OXT':['C'],'ND1':['CG'],'NE2':['CE1']},
              'HIS':{'O':['C'],'OXT':['C'],'NE2':['CE1'],'ND1':['CG']},
              'HIP':{'O':['C'],'OXT':['C'],'NE2':['CE1'],'ND1':['CG']},
              'TRP':{'O':['C'],'OXT':['C']},
              'PHE':{'O':['C'],'OXT':['C']},
              'TYR':{'O':['C'],'OXT':['C'],'OH':['CZ']},
              'GLU':{'O':['C'],'OXT':['C'],'OE1':['CD'],'OE2':['CD']},
              'ASP':{'O':['C'],'OXT':['C'],'OD1':['CG'],'OD2':['CG']},
              'LYS':{'O':['C'],'OXT':['C']},
              'LYN':{'O':['C'],'OXT':['C'],'NZ':['CE']},
              'CYS':{'O':['C'],'OXT':['C'],'SG':['CB']},
              'CYM':{'O':['C'],'OXT':['C'],'SG':['CB']},
              'CYX':{'O':['C'],'OXT':['C']},
              'MET':{'O':['C'],'OXT':['C']},
              'PRO':{'O':['C'],'OXT':['C']},
              'ASH':{'O':['C'],'OXT':['C'],'OD1':['CG']},
              'GLH':{'O':['C'],'OXT':['C'],'OE1':['CD']}
              }

    D_mappings={'ALA':{'N':['H']},
              'GLY':{'N':['H']},
              'SER':{'N':['H'],'OG':['HG']},
              'THR':{'N':['H'],'OG1':['HG1']},
              'LEU':{'N':['H'],},
              'ILE':{'N':['H'],},
              'VAL':{'N':['H'],},
              'ASN':{'N':['H'],'ND2':['HD21','HD22']},
              'GLN':{'N':['H'],'NE2':['HE21','HE22']},
              'ARG':{'N':['H'],'NE':['HE'],'NH1':['HH11','HH12'],'NH2':['HH21','HH22']},
              'HID':{'N':['H'],'ND1':['HD1']},
              'HIE':{'N':['H'],'NE2':['HE2']},
              'HIS':{'N':['H'],'ND1':['HD1'],'NE2':['HE2']},
              'HIP':{'N':['H'],'ND1':['HD1'],'NE2':['HE2']},
              'TRP':{'N':['H'],'NE1':['HE1']},
              'PHE':{'N':['H']},
              'TYR':{'N':['H'],'OH':['HH']},
              'GLU':{'N':['H']},
              'ASP':{'N':['H']},
              'LYS':{'N':['H'],'NZ':['HZ1','HZ2','HZ3']},
              'LYN':{'N':['H'],'NZ':['HZ2','HZ3']},
              'CYS':{'N':['H'],'SG':['HG']},
              'CYM':{'N':['H']},
              'CYX':{'N':['H']},
              'MET':{'N':['H']},
              'ASH':{'N':['H'],'OD2':['HD2']},
              'GLH':{'N':['H'],'OE2':['HE2']}
              }


    if prot_dict[p_coord[0]][p_coord[1]][-1] in ['P','D'] and p_coord[1] not in ['NZ']:
        p1=np.array([float(p_coord[2]),float(p_coord[3]),float(p_coord[4])])
        p2=np.array([float(b_coord[2]),float(b_coord[3]),float(b_coord[4])])
        if GetAtomDistance(p1,p2) <= 4.2:
            ang_list=[]
            for h in D_mappings[p_coord[0]][p_coord[1]]:
                h_coord=parsed_dict[p_coord[-1]][h]
                c=np.array([float(h_coord[2]),float(h_coord[3]),float(h_coord[4])])
                angle=get_angle(p1,p2,c)
                ang_list.append(angle)
            if max(ang_list) >= 90.0:
                return True,max(ang_list)
            else:
                return False,max(ang_list)
        else:
            return False, GetAtomDistance(p1,p2) 
    elif p_coord[0] in ['LYS','LYN'] and p_coord[1]=='NZ':
        c=np.array([float(p_coord[2]),float(p_coord[3]),float(p_coord[4])])
        p2=np.array([float(b_coord[2]),float(b_coord[3]),float(b_coord[4])])
        if GetAtomDistance(c,p2) <= 4.2:
            ce_coord=parsed_dict[p_coord[-1]]['CE']
            p1=np.array([float(ce_coord[2]),float(ce_coord[3]),float(ce_coord[4])])
            angle=get_angle(p1,p2,c)
            if angle >= 90.0:
                return True,angle
            else:
                return False,angle
        else:
            return False, GetAtomDistance(c,p2)
    elif prot_dict[p_coord[0]][p_coord[1]][-1] in ['A','N']:
        c=np.array([float(p_coord[2]),float(p_coord[3]),float(p_coord[4])])
        p2=np.array([float(b_coord[2]),float(b_coord[3]),float(b_coord[4])])
        if GetAtomDistance(c,p2) <= 4.2:
            c_coord=parsed_dict[p_coord[-1]][A_mappings[p_coord[0]][p_coord[1]][0]]
            p1=np.array([float(c_coord[2]),float(c_coord[3]),float(c_coord[4])])
            angle=get_angle(p1,p2,c)
            if angle >= 90.0:
                return True,angle
            else:
                return False,angle
        else:
            return False, GetAtomDistance(c,p2)
    elif prot_dict[p_coord[0]][p_coord[1]][-1] == 'DA' and p_coord[0] not in ['HIE','HID','HIP','HIS']:
        p1=np.array([float(p_coord[2]),float(p_coord[3]),float(p_coord[4])])
        p2=np.array([float(b_coord[2]),float(b_coord[3]),float(b_coord[4])])
        if GetAtomDistance(p1,p2) <= 4.2:
            p1a_coord=parsed_dict[p_coord[-1]][A_mappings[p_coord[0]][p_coord[1]][0]]
            p1d_coord=parsed_dict[p_coord[-1]][D_mappings[p_coord[0]][p_coord[1]][0]]
            p1a=np.array([float(p1a_coord[2]),float(p1a_coord[3]),float(p1a_coord[4])])
            p1d=np.array([float(p1d_coord[2]),float(p1d_coord[3]),float(p1d_coord[4])])
            angle_a=get_angle(p1a,p2,p1)
            angle_d=get_angle(p1,p2,p1d)
            if angle_a >= 90.0 or angle_d >= 90.0:
                return True,angle_a,angle_d
            else:
                return False,angle_a, angle_d
        else:
            return False, GetAtomDistance(p1,p2)
    elif prot_dict[p_coord[0]][p_coord[1]][-1] == 'DA' and p_coord[0] in ['HIE','HID','HIP','HIS']:
        c=np.array([float(p_coord[2]),float(p_coord[3]),float(p_coord[4])])
        p2=np.array([float(b_coord[2]),float(b_coord[3]),float(b_coord[4])])
        if GetAtomDistance(c,p2) <= 4.2:
            ce_coord=parsed_dict[p_coord[-1]][A_mappings[p_coord[0]][p_coord[1]][0]]
            p1=np.array([float(ce_coord[2]),float(ce_coord[3]),float(ce_coord[4])])
            angle=get_angle(p1,p2,c)
            if angle >= 85.0:
                return True,angle
            else:
                return False,angle
        else:
            return False, GetAtomDistance(c,p2) 
    elif prot_dict[p_coord[0]][p_coord[1]][-1] in ['AR','HP']:
        return False

def checking_free_polar_sites_witH(prot_coords,parsed_dict):
#hbd_list=['H','HG','HG1','HD21','HD22','HE21','HE22','HE','HH11','HH12','HH21','HH22','HD1','HE2','HH','HE1','HZ1','HZ2','HZ3','HE2','HD2']
#hba_list=['OG','OH','O','OD1','OD2','OE1','OE2','SG','SD','ND1']
    hba_list=['O','OP']
    hbd_list=['H','HP']

    for ix,i in enumerate(prot_coords):
        if i[1] in hbd_list:
            i_coords=[float(i[2]),float(i[3]),float(i[4])]
            for j in parsed_dict.keys():
                for k in parsed_dict[j].keys():
                    if k in hba_list:
                        k_coords=[float(parsed_dict[j][k][2]),float(parsed_dict[j][k][3]),float(parsed_dict[j][k][4])]
                        if GetAtomDistance(k_coords,i_coords) < 2.00 and parsed_dict[j][k][-1] != i[-1]:
#                            for lx,ll in enumerate(prot_coords):
#                                if ll[-2]==i[-2] and ll[1]=='N':
#                                    prot_coords[lx][1]=='NP'
#                                    print ll
                            prot_coords[ix][1]='HP'
#                            print i,parsed_dict[j][k],GetAtomDistance(k_coords,i_coords)
        elif i[1] in hba_list:
            i_coords=[float(i[2]),float(i[3]),float(i[4])]
            for j in parsed_dict.keys():
                for k in parsed_dict[j].keys():
                    if k in hbd_list:
                        k_coords=[float(parsed_dict[j][k][2]),float(parsed_dict[j][k][3]),float(parsed_dict[j][k][4])]
                        if GetAtomDistance(k_coords,i_coords) < 2.00 and parsed_dict[j][k][-1] != i[-1]:
                            prot_coords[ix][1]='OP'
#                            print i,parsed_dict[j][k],GetAtomDistance(k_coords,i_coords)

    for ix,i in enumerate(prot_coords):
        if i[1]=='HP':
            for lx,ll in enumerate(prot_coords):
                if ll[-1]==i[-1] and ll[1]=='N':
                    prot_coords[lx][1]='NP'
#                    print ll
    return prot_coords



def get_sasa_dict(prot_asa,pdb_directory):
    prot_apo_list=[]
    file_list=sorted(glob.glob(pdb_directory.strip('/')+'/*.pdb'))
    file_list.pop(-1)
    for i,file in enumerate(file_list):
        if file.endswith("beta.pdb"):
            file_list.pop(i)
    for filename in file_list:
        with open(filename,'r') as f:
            templines=f.readlines()
        for i in templines:
            if i[17:20]=='AAC':
                break
            else:
                prot_apo_list.append([i[22:26].strip(),i[12:16].strip()])
    prot_apo_list.sort()
    prot_apo_list=list(prot_apo_list for prot_apo_list,_ in itertools.groupby(prot_apo_list))

    with open(glob.glob(prot_asa)[0],'r') as f:
        prot_apo_lines=f.readlines()

    prot_apo_surf=defaultdict(dict)
    for i in prot_apo_lines:
        res_id=i[22:26].strip()
        atom_name=i[12:16].strip()
        prot_apo_surf[res_id][atom_name]=i

    sasa_dict=defaultdict(dict)
    for i in prot_apo_list:
        apo_surf=float(prot_apo_surf[i[0]][i[1]][56:62:].strip())
        sasa_dict[prot_apo_surf[i[0]][i[1]][22:26].strip()][prot_apo_surf[i[0]][i[1]][12:16].strip()]=[prot_apo_surf[i[0]][i[1]][17:20],prot_apo_surf[i[0]][i[1]][12:16].strip(),apo_surf,prot_apo_surf[i[0]][i[1]][22:26].strip()]

    return sasa_dict

def get_beta_prot_sa(prot_pdb,beta_pdb):
    with open(prot_pdb,'r') as f:
        prot_lines=f.readlines()
    with open(beta_pdb,'r') as f:
        beta_lines=f.readlines()
    prot_lines.append("TER\n")
    for k in beta_lines:
        if k[13:16]=='BAO':
            prot_lines.append(k)
    prot_lines.append("\n")
    prot_lines.append("TER\n")
    with open("temp_beta.pdb","w") as g:
        g.writelines("%s" % item for item in prot_lines)
    os.system("/Users/jkatigba/bin/Naccess/naccess temp_beta.pdb")
    with open("temp_beta.asa",'r') as f:
        prot_beta=f.readlines()
    beta_surf=[]
    prot_beta_list=[]
    for i in prot_beta:
        if i[17:20]=='BAC':
            beta_surf.append([i[17:20],i[13:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),float(i[56:62].strip())])
        else:
            res_id=i[22:26].strip()
            atom_name=i[13:16].strip()
            prot_beta_list.append(i)
    os.remove("temp_beta.pdb")
    os.remove("temp_beta.rsa")
    os.remove("temp_beta.log")
    os.remove("temp_beta.asa")
    return beta_surf, prot_beta_list

def get_lig_metal_coordinates_from_ac(pdb_file):
    metal_list=['Na','Mg','K','Ca','V','Mn','Fe','Co','Ni','Cu','Zn','Mo','W','Cd','Hg']
    with open(pdb_file,'r') as f:
        templines=f.readlines()
    lig_temp=[]
    metal_temp=[]
    for k in templines:
        if k[76:78] in metal_list:
            metal_temp.append(k)
        else:
            lig_temp.append(k)
    with open('temp_a.pdb','w') as h:
        h.writelines(lig_temp)
    os.system("obabel temp_a.pdb -O temp.pdb -p 7")
    os.system("~/amber14/bin/antechamber -fi pdb -i temp.pdb -o lig.ac -fo ac")
    with open('lig.ac','r') as f:
        templines=f.readlines()
    Coordinates=[]
    templines=templines[2:]
    for i in templines:
        if i[12:16].strip()=='H':
            break
        elif i[0:4]=='BOND':
            break
        else:
            Coordinates.append(['LIG',i[72:].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
    for j in metal_temp:
        Coordinates.append(['LIG','m',float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])

    os.system('rm ANTECHAMBER* ATOMTYPE* lig.ac temp.pdb temp_a.pdb')
    return Coordinates

## writes the lig coords and splits the metal ions from the lig
def write_lig_coords_pdb(pdb_file):
    metal_list=['Na','Mg','K','Ca','V','Mn','Fe','Co','Ni','Cu','Zn','Mo','W','Cd','Hg']
    with open(pdb_file,'r') as f:
        templines=f.readlines()
    lig_temp=[]
    metal_temp=[]
    for k in templines:
        if k[76:78] in metal_list:
            metal_temp.append(k)
        else:
            lig_temp.append(k)
    if metal_temp:
        with open('metal.pdb','w') as h:
            h.writelines(metal_temp)
    with open('lig.pdb','w') as f:
        f.writelines(lig_temp)


def get_coordinates_pdb_withID(pdb_file):
    Coordinates=[]
    with open(pdb_file,'r') as f:
        templines=f.readlines()
    for i in templines:
        Coordinates.append([i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip()),i[22:26].strip()])
    return Coordinates

def get_coordinates_pdb(pdb_file):
    Coordinates=[]
    with open(pdb_file,'r') as f:
        templines=f.readlines()
    for i in templines:
        Coordinates.append([i[17:20].strip(),i[12:16].strip(),float(i[30:38].strip()),float(i[38:46].strip()),float(i[46:54].strip())])
    return Coordinates

def get_parsed_coord_list(prot_coords,parsed_dict):
    res_list=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','HID','HIE','HIP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','ASH','GLH','CYM','CYX','LYN']
    pdb_coords=[]
    pdb_data=[]
    for i in prot_coords:
        if i[0] not in res_list:
            pdb_coords.append([float(i[2]),float(i[2]),float(i[3])])
            pdb_data.append([i[0],i[1],i[-1]])
    for rd in parsed_dict.keys():
        for ad in parsed_dict[rd].keys():
            if ad[0]!='H':
                pdb_coords.append([float(parsed_dict[rd][ad][2]),float(parsed_dict[rd][ad][3]),float(parsed_dict[rd][ad][4])])
                pdb_data.append([parsed_dict[rd][ad][0],parsed_dict[rd][ad][1],parsed_dict[rd][ad][-1]])
    return np.array(pdb_coords),np.array(pdb_data)

def transform_grid(trans_coords,untrans_coords,grid_coords):
#def transform_grid(trans_coords,untrans_coords,grid_coords,out_file):
    A=[]
    for i in untrans_coords:
        A.append([float(i[2]),float(i[3]),float(i[4])])
    B=[]
    for j in trans_coords:
        B.append([float(j[2]),float(j[3]),float(j[4])])
    A=mat(A)
    B=mat(B)
    assert len(A) == len(B)

    N = A.shape[0]; # total points

    centroid_A = mean(A, axis=0)
    centroid_B = mean(B, axis=0)

    # centre the points
    AA = A - tile(centroid_A, (N, 1))
    BB = B - tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = transpose(AA) * BB

    U, S, Vt = linalg.svd(H)

    R = Vt.T * U.T

    # special reflection case
    if linalg.det(R) < 0:
       print("Reflection detected")
       Vt[2,:] *= -1
       R = Vt.T * U.T

    t = -R*centroid_A.T + centroid_B.T
    C=[]
    for k in grid_coords:
        C.append([float(k[2]),float(k[3]),float(k[4])])
    D=(R*mat(C).T)+tile(t,(1,len(C)))
    D=D.T
    x=len(D)
    transform=[]
    for lx,l in enumerate(grid_coords):
        transform.append([l[0],l[1],round(D[lx,0],3),round(D[lx,1],3),round(D[lx,2],3)])
    return transform


def fill_p_verts(p_verts,cluster):
    """
    Function: fill an empty list of lists (length equal to the number of clusters) with vertice indices according to the clustering (becomes centroids in super clustering)
    Parameters
    p_verts: empty list of lists
    cluster: cluster indices 
    """

    for i,clust_idx in enumerate(cluster):
        p_idx = clust_idx - 1
        p_verts[p_idx].append(i)
    return

def get_pharm_clust_coord(grid_coords,clust_dist):
    g_coords = []
    for g in grid_coords:
        g_coords.append([g[2],g[3],g[4]])
    zmat = hier.linkage(g_coords, method='average')
    #zmat = hier.linkage(g_coords, method='complete')
    cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
    #cluster = hier.fcluster(zmat,3.5,criterion='distance')

    clusters = [ [] for i in range(max(cluster)) ]
    fill_p_verts(clusters,cluster)

    frag_centroids = []

    for count,i in enumerate(clusters):
        frag_coords = []
        for idx in i:
            frag_coords.append(g_coords[idx])
        frag_centroid = np.average(frag_coords,axis=0)

        frag_centroids.append(frag_centroid)

    type_name=grid_coords[0][1]
    new_coords=[]
    for jx,jj in enumerate(frag_centroids):
        new_coords.append(['GRD',type_name,round(jj[0],3),round(jj[1],3),round(jj[2],3)])
#    for jx,jj in enumerate(frag_centroids):
#        new_coords.append('HETATM'+str(jx).rjust(5)+'  H           0    '+str(round(jj[0],3)).rjust(8)+str(round(jj[1],3)).rjust(8)+str(round(jj[2],3)).rjust(8)+'  0.00  0.00           C\n')

    return new_coords



#def get_pharm_clust_coord_screen(grid_coords,clust_dist,mem_thresh):
#    g_coords = []
#    for g in grid_coords:
#        g_coords.append([g[2],g[3],g[4]])
#    zmat = hier.linkage(g_coords, method='average')
#    #zmat = hier.linkage(g_coords, method='complete')
#    cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
#    #cluster = hier.fcluster(zmat,3.5,criterion='distance')
#
#    clusters = [ [] for i in range(max(cluster)) ]
#    fill_p_verts(clusters,cluster)
#
#    frag_centroids = []
#
#    for count,i in enumerate(clusters):
#        frag_coords = []
#        if len(i) > mem_thresh:
#            for idx in i:
#                frag_coords.append(g_coords[idx])
#        frag_centroid = np.average(frag_coords,axis=0)
#
#        frag_centroids.append(frag_centroid)
#
#    type_name=grid_coords[0][1]
#    new_coords=[]
#    for jx,jj in enumerate(frag_centroids):
#        new_coords.append(['GRD',type_name,round(jj[0],3),round(jj[1],3),round(jj[2],3)])
#    for jx,jj in enumerate(frag_centroids):
#        new_coords.append('HETATM'+str(jx).rjust(5)+'  H           0    '+str(round(jj[0],3)).rjust(8)+str(round(jj[1],3)).rjust(8)+str(round(jj[2],3)).rjust(8)+'  0.00  0.00           C\n')

    return new_coords

if __name__=="__main__":
    print(res)

