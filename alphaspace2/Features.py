import math
import os
from collections import defaultdict

import numpy as np
import scipy
from scipy import stats
from scipy.spatial import Voronoi
from scipy.spatial.distance import cdist
from scipy.special import cbrt

from alphaspace2.functions import getSASA

#### Default mass of beta atoms is 12.0 (Same as carbon)
_c_mass = 12.0
_c2_mass = 12.0 * 12.0


####

def _get_types_dict(filename):
    types_dict = defaultdict(dict)
    fp = open(filename, 'r')
    lines = fp.readlines()
    fp.close()
    for line in lines:
        line = line.split()
        types_dict[line[2]][line[0]] = line

    return types_dict


atom_types = ['P', 'N', 'D', 'A', 'DA', 'AR', 'HP', 'PL']
atom_types_match = {'P': 'Positive_OASA', 'N': 'Negative_OASA', 'D': 'H_bond_Donor_OASA', 'A': 'H_bond_Acceptor_OASA',
                    'DA': 'H_bond_Doneptor_OASA', 'AR': 'Aromatic_OASA', 'HP': 'Hydrophobic_OASA', 'PL': 'Polar_OASA',
                    'NL': 'Null_type_OASA'}
this_dir, this_filename = os.path.split(__file__)
prot_dict = _get_types_dict(this_dir + '/data/amber_types_fine.dat')
res_names = prot_dict.keys()


def _get_grid_volume(coords_array, dd=0.5, radii=1.6):
    """

    :param coords_array: n*3
    :param dd:
    :param radii:
    :return: float
    """
    beta_temp_coord = coords_array - np.mean(coords_array, axis=0)
    if len(beta_temp_coord) > 1:
        x_min, y_min, z_min = np.min(beta_temp_coord, axis=0)
        x_max, y_max, z_max = np.max(beta_temp_coord, axis=0)
    else:
        x_min, y_min, z_min = beta_temp_coord[0][0], beta_temp_coord[0][1], beta_temp_coord[0][2]
        x_max, y_max, z_max = beta_temp_coord[0][0], beta_temp_coord[0][1], beta_temp_coord[0][2]

    x_grid = np.arange(x_min - radii, x_max + radii + 0.1, dd)
    y_grid = np.arange(y_min - radii, y_max + radii + 0.1, dd)
    z_grid = np.arange(z_min - radii, z_max + radii + 0.1, dd)

    grid_coords = np.vstack(np.meshgrid(x_grid, y_grid, z_grid)).reshape(3, -1).T

    dist = scipy.spatial.distance.cdist(beta_temp_coord, grid_coords)

    distance_compare = (dist[0] <= radii)
    for d in range(1, len(dist)):
        distance_compare = np.logical_or(distance_compare, (dist[d] <= radii))

    return pow(dd, 3) * distance_compare.sum()


def _get_overlap_volume(coords_array_1, coords_array_2, dd=0.5, radii=1.6):
    """

    :param coords_array_1:
    :param coords_array_2:
    :param dd:
    :param radii:
    :return: float
    """
    cent_point = np.mean(coords_array_1)
    coords_array_1 = coords_array_1 - cent_point
    coords_array_2 = coords_array_2 - cent_point

    if len(coords_array_1) > 1:
        x_min_1, y_min_1, z_min_1 = np.min(coords_array_1, axis=0)
        x_max_1, y_max_1, z_max_1 = np.max(coords_array_1, axis=0)
    else:
        x_min_1, y_min_1, z_min_1 = coords_array_1[0][0], coords_array_1[0][1], coords_array_1[0][2]
        x_max_1, y_max_1, z_max_1 = coords_array_1[0][0], coords_array_1[0][1], coords_array_1[0][2]

    if len(coords_array_2) > 1:
        x_min_2, y_min_2, z_min_2 = np.min(coords_array_2, axis=0)
        x_max_2, y_max_2, z_max_2 = np.max(coords_array_2, axis=0)
    else:
        x_min_2, y_min_2, z_min_2 = coords_array_2[0][0], coords_array_2[0][1], coords_array_2[0][2]
        x_max_2, y_max_2, z_max_2 = coords_array_2[0][0], coords_array_2[0][1], coords_array_2[0][2]

    x_min, y_min, z_min = min(x_min_1, x_min_2), min(y_min_1, y_min_2), min(z_min_1, z_min_2)
    x_max, y_max, z_max = max(x_max_1, x_max_2), max(y_max_1, y_max_2), max(z_max_1, z_max_2)

    x_grid = np.arange(x_min - radii, x_max + radii + 0.1, dd)
    y_grid = np.arange(y_min - radii, y_max + radii + 0.1, dd)
    z_grid = np.arange(z_min - radii, z_max + radii + 0.1, dd)
    grid_coords = np.vstack(np.meshgrid(x_grid, y_grid, z_grid)).reshape(3, -1).T
    dist_1 = scipy.spatial.distance.cdist(coords_array_1, grid_coords)
    dist_2 = scipy.spatial.distance.cdist(coords_array_2, grid_coords)
    distance_compare_1 = (dist_1[0] <= 1.6)
    distance_compare_2 = (dist_2[0] <= 1.6)
    for d in range(1, len(dist_1)):
        distance_compare_1 = np.logical_or(distance_compare_1, (dist_1[d] <= radii))
    for d in range(1, len(dist_2)):
        distance_compare_2 = np.logical_or(distance_compare_2, (dist_2[d] <= radii))
    return pow(dd, 3) * np.logical_and(distance_compare_1, distance_compare_2).sum()


def _soergel(i, j):
    """

    :param i: np array
    :param j: np array
    :return: float
    """
    ### Distance calculation (Soergel Distance)
    if len(i) is not len(j):
        raise('lengths not equal')
    else:
        ntemp = 0
        dtemp = 0
        for k in range(len(i)):
            ntemp = ntemp + abs(i[k] - j[k])
            dtemp = dtemp + max(i[k], j[k])
        score = float(ntemp) / float(dtemp)
    return score


def _get_pharmacophore_fingerprint(prot_md, beta_atoms):
    """

    :param prot_md:
    :param beta_atoms:
    :return: dict
    """
    s1 = getSASA(prot_md)
    s2 = getSASA(prot_md, cover_atom_coords=beta_atoms / 10)
    diff_bool = s1 - s2 > 0
    d_asa = (s1 - s2)[diff_bool] * 100
    topology = prot_md.topology
    table, bonds = topology.to_dataframe()
    occluded_asa = []
    for top, asa in zip(table.values[diff_bool], d_asa):
        if top[4] in res_names:
            if top[1] in prot_dict[top[4]].keys():
                occluded_asa.append(tuple([prot_dict[top[4]][top[1]][-1], asa]))
            else:
                occluded_asa.append(tuple(['NL', asa]))
        else:
            occluded_asa.append(tuple(['NL', asa]))
    pharma_fingerprint_dict = {'Total_OASA': 0, 'Positive_OASA': 0, 'Negative_OASA': 0, 'H_bond_Donor_OASA': 0,
                      'H_bond_Acceptor_OASA': 0,
                      'H_bond_Doneptor_OASA': 0, 'Aromatic_OASA': 0, 'Hydrophobic_OASA': 0, 'Polar_OASA': 0,
                      'Null_type_OASA': 0}
    for dd in occluded_asa:
        pharma_fingerprint_dict['Total_OASA'] += dd[1]
        pharma_fingerprint_dict[atom_types_match[dd[0]]] += dd[1]
    return pharma_fingerprint_dict


def _Get_USR_alpha_beta(coords):
    """

    :param coords:
    :return: list
    """
    def _ctd_alpha(coords):
        """

        :param coords: n * 3
        :return:
        """
        ctd = np.array([np.mean(coords, axis=0)])
        ctd_distance = cdist(ctd, coords)[0]
        ctd_alpha = {}
        ctd_alpha['ctd_1'] = round(np.mean(ctd_distance), 3)
        ctd_alpha['ctd_2'] = round(np.std(ctd_distance), 3)
        ctd_alpha['ctd_3'] = round(cbrt(stats.skew(ctd_distance)), 3)
        return ctd_alpha

    def _cst_alpha(coords):
        ctd = np.array([np.mean(coords, axis=0)])
        ctd_distance = cdist(ctd, coords)[0]

        cst_index = np.argmin(ctd_distance)
        cst = np.array([coords[cst_index]])
        cst_distance = cdist(cst, coords)[0]

        cst_alpha = {}
        cst_alpha['cst_1'] = round(np.mean(cst_distance), 3)
        cst_alpha['cst_2'] = round(np.std(cst_distance), 3)
        cst_alpha['cst_3'] = round(cbrt(stats.skew(cst_distance)), 3)
        return cst_alpha

    def _fct_alpha(coords):
        ctd = np.array([np.mean(coords, axis=0)])
        ctd_distance = cdist(ctd, coords)[0]

        fct_index = np.argmax(ctd_distance)
        fct = np.array([coords[fct_index]])
        fct_distance = cdist(fct, coords)[0]

        fct_alpha = {}
        fct_alpha['fct_1'] = round(np.mean(fct_distance), 3)
        fct_alpha['fct_2'] = round(np.std(fct_distance), 3)
        fct_alpha['fct_3'] = round(cbrt(stats.skew(fct_distance)), 3)
        return fct_alpha

    def _ftf_alpha(coords):
        ctd = np.array([np.mean(coords, axis=0)])
        ctd_distance = cdist(ctd, coords)[0]

        fct_index = np.argmax(ctd_distance)
        fct = np.array([coords[fct_index]])
        fct_distance = cdist(fct, coords)[0]

        ftf_index = np.argmax(fct_distance)
        ftf = np.array([coords[ftf_index]])
        ftf_distance = cdist(ftf, coords)[0]

        ftf_alpha = {}
        ftf_alpha['ftf_1'] = round(np.mean(ftf_distance), 3)
        ftf_alpha['ftf_2'] = round(np.std(ftf_distance), 3)
        ftf_alpha['ftf_3'] = round(cbrt(stats.skew(ftf_distance)), 3)
        return ftf_alpha

    result = {}
    result.update(_ctd_alpha(coords))
    result.update(_cst_alpha(coords))
    result.update(_fct_alpha(coords))
    result.update(_ftf_alpha(coords))
    return np.array(list(result.values()))


def _calculate_3DWiener(DistanceMatrix):
    """
    :param DistanceMatrix:
    :return: float
    """
    return np.round((scipy.sum(DistanceMatrix) / 2.0) / 1000, 3)


def _calculate_Petitjean3DIndex(DistanceMatrix):
    temp1 = scipy.amax(DistanceMatrix, axis=0)
    temp2 = round(max(temp1) / min(temp1) - 1.0, 3)
    if math.isnan(temp2):
        return 0
    else:
        return temp2


def _calculate_GeometricDiameter(DistanceMatrix):
    temp1 = scipy.amax(DistanceMatrix, axis=0)
    return max(temp1)


def _calculate_Harary3D(DistanceMatrix):
    nAT = len(DistanceMatrix)
    res = 0.0
    for i in range(nAT - 1):
        for j in range(i + 1, nAT):
            if DistanceMatrix[i, j] == 0:
                cds = 0.0
            else:
                cds = 1. / DistanceMatrix[i, j]
            res = res + cds
    return res/10


def _calculate_AverageGeometricDistanceDegree(DistanceMatrix):
    nAT = len(DistanceMatrix)
    res = sum(sum(DistanceMatrix)) / nAT
    return res


def _calculate_Gravitational3D1(Coordinates):
    nAT = len(Coordinates)
    result = 0.0
    for i in range(nAT - 1):
        for j in range(i + 1, nAT):
            dis = np.linalg.norm(Coordinates[i] - Coordinates[j])
            result = result + _c_mass / scipy.power(dis, p=2)
    return result / 100.


def _calculate_RadiusofGyration(coords):
    mass = 12.0
    masscenter = np.mean(coords, axis=0)
    result = 0.0
    for i in range(len(coords)):
        dis = np.linalg.norm(np.abs(coords[i] - masscenter))
        result = result + _c_mass * scipy.power(dis, p=2)
    return scipy.sqrt(result / mass)


def _calculate_SPAN(coords):
    masscenter = np.array([np.mean(coords, axis=0)])
    res = cdist(coords, masscenter).flatten()
    return max(res)


def _GetInertiaMatrix(Coordinates):
    masscenter = np.mean(Coordinates, axis=0)
    Coordinates = Coordinates - masscenter
    nAT = len(Coordinates)
    InertiaMatrix = scipy.zeros((3, 3))
    res11, res22, res33, res12, res23, res13 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    for i in range(nAT):
        res11 = res11 + _c_mass * (np.power(Coordinates[i][1], 2) + np.power(Coordinates[i][2], 2))
        res22 = res22 + _c_mass * (np.power(Coordinates[i][0], 2) + np.power(Coordinates[i][2], 2))
        res33 = res33 + _c_mass * (np.power(Coordinates[i][0], 2) + np.power(Coordinates[i][1], 2))
        res12 = res12 + _c_mass * (Coordinates[i][0] * Coordinates[i][1])
        res13 = res13 + _c_mass * (Coordinates[i][0] * Coordinates[i][2])
        res23 = res23 + _c_mass * (Coordinates[i][1] * Coordinates[i][2])
    InertiaMatrix[0, 0], InertiaMatrix[1, 1], InertiaMatrix[2, 2] = res11, res22, res33
    InertiaMatrix[0, 1], InertiaMatrix[1, 0] = -res12, -res12
    InertiaMatrix[0, 2], InertiaMatrix[2, 0] = -res13, -res13
    InertiaMatrix[1, 2], InertiaMatrix[2, 1] = -res23, -res23
    return InertiaMatrix


def _calculate_MolecularEccentricity(Coordinates):
    InertiaMatrix = _GetInertiaMatrix(Coordinates)
    u, s, v = scipy.linalg.svd(InertiaMatrix)
    res1 = s[0]
    res3 = s[2]

    res = math.pow(res1 * res1 - res3 * res3, 1. / 2) / res1
    if math.isnan(res):
        return 0
    else:
        return res


def _calculate_PrincipalMomentofInertia(Coordinates):
    InertiaMatrix = _GetInertiaMatrix(Coordinates)
    u, s, v = scipy.linalg.svd(InertiaMatrix)
    res = {}
    res['IA'] = round(s[2], 3)
    res['IB'] = round(s[1], 3)
    res['IC'] = round(s[0], 3)
    return res


def _calculate_RatioPMI(Coordinates):
    temp = _calculate_PrincipalMomentofInertia(Coordinates)
    res = {}
    if temp['IB'] == 0 or temp['IC'] == 0:
        res['IA/B'] = 0
        res['IA/C'] = 0
        res['IB/C'] = 0
    else:
        res['IA/B'] = round(temp['IA'] / temp['IB'], 3)
        res['IA/C'] = round(temp['IA'] / temp['IC'], 3)
        res['IB/C'] = round(temp['IB'] / temp['IC'], 3)
    return res


def _calculate_NormalizedRatioPMI(Coordinates):
    pmi_tuples = [(b, a) for a, b in _calculate_PrincipalMomentofInertia(Coordinates).items()]
    pmi_tuples = sorted(pmi_tuples)
    res = {}
    if pmi_tuples[1][0] == 0 or pmi_tuples[2][0] == 0:
        res['NPR1'] = 0
        res['NPR2'] = 0
    else:
        res['NPR1'] = round(pmi_tuples[0][0] / pmi_tuples[2][0], 3)
        res['NPR2'] = round(pmi_tuples[1][0] / pmi_tuples[2][0], 3)
    return res
