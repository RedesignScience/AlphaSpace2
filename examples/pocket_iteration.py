from alphaspace import *
import sys

"""load processed .as file"""
universe = load_universe(sys.argv[1])

"""Run alphaspace with 4 cpu cores"""
universe.run_alphaspace_mp()

"""Iterate over pockets in a particular snapshot"""
snapshot_idx = 0
pockets = list(universe.pockets(snapshot_idx))
pockets.sort(key=lambda p: p.space, reverse=True)

for pocket in pockets:
    print("_Pocket Space {}".format(pocket.space))
    for beta in pocket.betas:
        print("_Beta Space {:.2f} at {}".format(beta.space, beta.centroid))
