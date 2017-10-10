import numpy as np
from scipy.spatial.distance import cdist
from AS.AS_Config import AS_Config


def getTetrahedronVolume(coord_list):
    """
    generate the volume
    :param coord_list: list N*4*3
    :return: float
    """
    coord_matrix = np.concatenate((np.array(coord_list),np.ones((4,1))),axis=1)

    return np.abs(np.linalg.det(coord_matrix) / 6)


def checkContact(coord_list_1,coord_list_2,threshold=None):
    """
    check which one in the coordinate list is in contact with the second coordinate
    :param coord_list_1: list of array N*3
    :param coord_list_2: list of array M*3
    :param threshold: float
    :return: array
    """
    if threshold is None:
        threshold = AS_Config().contact_threshold
    print(coord_list_1.shape,coord_list_2.shape)
    distance_matrix = cdist(coord_list_1,coord_list_2)
    return np.where(distance_matrix < threshold)

# def ifContact()



def getGridVolume(coord_list,threshold=1.6,resolution= 0.05):

    """
    calculate the volume of a point set using grid point approximation
    :param coord_list: array N * 3
    :param threshold: float
    :param resolution: float
    :return: float
    """
    max_coord = np.max(coord_list,axis=0)
    min_coord = np.min(coord_list,axis=0)
    coord_range = np.array([min_coord,max_coord]).transpose().tolist()
    x,y,z = [np.arange(start=ax[0]-threshold,stop=ax[1]+threshold,step=resolution) for ax in coord_range]
    grid_coords = np.array(np.meshgrid(x,y,z)).transpose().reshape((-1,3))

    grid_count = len(checkContact(grid_coords,coord_list,threshold=threshold)[0])
    return grid_count*(resolution**3)


