import numpy as np
from scipy.spatial import distance


class dPocket(object):
    def __init__(self, sPocket=[]):
        self.sPocket = sPocket
        self.sPocket_centers = []

    def get_sPockets(self, idx=()):
        if idx == ():
            return self.sPocket
        else:
            return (self.sPocket[i] for i in idx)

    def add_sPockets(self, pocket_list):
        self.sPocket.extend(pocket_list)

    def add_sPocket(self, pocket):
        self.sPocket.append(pocket)

    def get_sPocket_centers(self):
        if len(self.sPocket_centers) == len(self.sPocket):
            return self.sPocket_centers
        else:
            self.sPocket_centers = np.array([pocket.vert_cntr() for pocket in self.sPocket])
            return self.get_sPocket_centers()

    def __sub__(self, other):
        d = np.average(distance.cdist(self.get_sPocket_centers(), other.get_sPocket_centers()))
        return d

    def get_scores(self):
        return [s.score for s in self.sPocket]

    def get_average_score(self):
        return self.get_sum_score()/len(self)

    def get_sum_score(self):
        return np.sum(self.get_scores())

    def merge(self, other):
        self.sPocket = [p1 for p1 in self.get_sPockets()]+[p2 for p2 in other.get_sPockets()]
        self.sPocket_centers = []

    def get_ss_idx(self):
        return set([i.ss_idx for i in self.sPocket])

    def __len__(self):
        return len(set([p.ss_idx for p in self.get_sPockets()]))

