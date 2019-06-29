import numpy as np
from .beta_usr import Get_USR_alpha_beta

## usr keys represent the different USR based features. USR features are the 1st, 2nd, and 3rd moments of the distribution of distances of atoms in the molecules to 4-specific points (the molecular centroid (ctd), the closest atom to ctd (cst), the farthest atom to ctd (fct), and the farthest atom to fct (ftf))

usr_keys=['cst_1','cst_2','cst_3','ctd_1','ctd_2','ctd_3','fct_1','fct_2','fct_3','ftf_1','ftf_2','ftf_3']

## example coordinates

alpha_atom_pdb="""
ATOM      0  AAO AAC H   0      56.383  24.574  27.371        3.98
ATOM      1  AAO AAC H   0      55.937  23.843  27.780        3.78 
ATOM      2  AAO AAC H   0      56.095  24.254  27.513        3.89
ATOM      3  AAO AAC H   0      56.152  24.143  27.517        3.92
ATOM      4  AAU AAC H   0      58.402  24.859  28.412        3.88
ATOM      5  AAU AAC H   0      58.641  24.731  29.857        3.25
ATOM      6  AAU AAC H   0      58.990  25.209  31.729        3.27
ATOM      7  AAO AAC H   0      56.778  24.463  27.280        4.17
ATOM      8  AAO AAC H   0      56.517  24.268  27.391        4.06
ATOM      9  AAU AAC H   0      54.846  23.854  28.513        3.25
ATOM     10  AAU AAC H   0      58.591  24.695  29.845        3.22
ATOM     11  AAO AAC H   0      56.439  23.512  27.878        3.68
ATOM     12  AAU AAC H   0      55.391  23.501  28.283        3.52
ATOM     13  AAU AAC H   0      54.893  23.368  28.672        3.37
ATOM     14  AAU AAC H   0      54.718  23.057  28.764        3.35
ATOM     15  AAU AAC H   0      54.708  23.079  28.776        3.34
ATOM     16  AAO AAC H   0      56.207  23.735  27.854        3.78
ATOM     17  AAO AAC H   0      56.442  23.520  27.873        3.68
ATOM     18  AAO AAC H   0      56.578  23.848  27.736        3.90
ATOM     19  AAO AAC H   0      57.955  24.737  26.794        4.39
ATOM     20  AAU AAC H   0      60.452  22.976  26.609        3.27
ATOM     21  AAO AAC H   0      57.395  24.812  25.632        3.85
ATOM     22  AAO AAC H   0      57.610  24.862  25.699        3.93
ATOM     23  AAO AAC H   0      57.292  24.616  26.304        4.00
ATOM     24  AAO AAC H   0      57.860  24.726  26.576        4.29
ATOM     25  AAO AAC H   0      58.205  24.265  26.523        3.96
ATOM     26  AAU AAC H   0      60.125  22.517  26.412        3.30
ATOM     27  AAO AAC H   0      60.060  22.873  26.781        3.52
ATOM     28  AAO AAC H   0      56.958  23.422  28.002        3.54
ATOM     29  AAO AAC H   0      59.882  22.887  27.324        3.51
ATOM     30  AAO AAC H   0      58.047  24.629  26.804        4.32
ATOM     31  AAO AAC H   0      57.169  24.046  27.757        3.97
ATOM     32  AAO AAC H   0      56.587  23.851  27.735        3.90
ATOM     33  AAO AAC H   0      56.487  23.558  27.851        3.71
ATOM     34  AAO AAC H   0      56.579  23.841  27.738        3.89
ATOM     35  AAO AAC H   0      57.847  24.652  27.383        4.42
ATOM     36  AAO AAC H   0      57.800  24.716  27.298        4.46
ATOM     37  AAO AAC H   0      57.479  24.613  27.281        4.36
ATOM     38  AAO AAC H   0      57.897  24.717  26.966        4.41
ATOM     39  AAU AAC H   0      60.310  22.823  27.243        3.31
ATOM     40  AAU AAC H   0      60.425  22.817  27.554        3.22
ATOM     41  AAO AAC H   0      60.047  22.825  27.005        3.51
ATOM     42  AAO AAC H   0      59.992  22.819  27.214        3.51
ATOM     43  AAU AAC H   0      60.125  22.517  26.412        3.30
"""


## from shape pdb coords and pass to Get_USR_alpha_beta
coords=np.array([[float(l[30:38].strip()),float(l[38:46].strip()),float(l[46:54].strip())] for l in alpha_atom_pdb.split('\n') if l[0:6] in ['ATOM  ','HETATM']])
bac_coords=[['XX','XXX',l[0],l[1],l[2]] for l in coords]
usr_dict=Get_USR_alpha_beta(bac_coords)
usr_array=[usr_dict[k] for k in usr_keys]

## the usr_array will be a feature vector with length 12. In my original implementation i would change the XX and XXX entries to represent different atom or residue types (eg ['GLU','OE1',59.992  22.819  27.214]) which can be used to calculate pharmacophoric-based USR features


## to get the similarity I use the Soergel distance

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
        score=float(ntemp)/float(dtemp)
    return(score)


## other distances from scipy.spatial can also be used, soergel is nice because the distance is between normalized between [0,1] and similarity is just 1-distance

similarity = 1 - soergel(usr_array_1,usr_array_2)




