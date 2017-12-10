"""
This function calculate the h-bond formed between the trajectories of two molecules
"""

import mdtraj
import numpy as np

mdtraj.baker_hubbard()


def baker_hubbard(traj1, traj2, freq=0.1, periodic=True):
    """Identify hydrogen bonds based on cutoffs for the Donor-H...Acceptor
    distance and angle.

    The criterion employed is :math:`\\theta > 120` and
    :math:`r_\\text{H...Acceptor} < 2.5 A`.

    When donor the donor is 'N' and the acceptor is 'O', this corresponds to
    the definition established in [1]_. The donors considered by this method
    are NH and OH, and the acceptors considered are O and N.

    Parameters
    ----------
    traj : md.Trajectory
        An mdtraj trajectory. It must contain topology information.
    freq : float, default=0.1
        Return only hydrogen bonds that occur in greater this fraction of the
        frames in the trajectory.
    exclude_water : bool, default=True
        Exclude solvent molecules from consideration
    periodic : bool, default=True
        Set to True to calculate displacements and angles across periodic box boundaries.
    sidechain_only : bool, default=False
        Set to True to only consider sidechain-sidechain interactions.

    Returns
    -------
    hbonds : np.array, shape=[n_hbonds, 3], dtype=int
        An array containing the indices atoms involved in each of the identified
        hydrogen bonds. Each row contains three integer indices, `(d_i, h_i,
        a_i)`, such that `d_i` is the index of the donor atom, `h_i` the index
        of the hydrogen atom, and `a_i` the index of the acceptor atom involved
        in a hydrogen bond which occurs (according to the definition above) in
        proportion greater than `freq` of the trajectory.

    Notes
    -----
    Each hydrogen bond is distinguished for the purpose of this function by the
    indices of the donor, hydrogen, and acceptor atoms. This means that, for
    example, when an ARG sidechain makes a hydrogen bond with its NH2 group,
    you might see what appear like double counting of the h-bonds, since the
    hydrogen bond formed via the H_1 and H_2 are counted separately, despite
    their "chemical indistinguishably"


    References
    ----------
    .. [1] Baker, E. N., and R. E. Hubbard. "Hydrogen bonding in globular
        proteins." Progress in Biophysics and Molecular Biology
        44.2 (1984): 97-179.
    """

    # Cutoff criteria: these could be exposed as function arguments, or
    # modified if there are better definitions than the this one based only
    # on distances and angles

    def _get_bond_triplets_between_structures(donor_top, acceptor_top):
        """
        Return the triplets of donor and acceptor in given order.
        :param donor_top: mdtraj.topology
        :param acceptor_top: mdtraj.topology
        :return:
        """

        def get_donors(e0, e1, top):
            # Find all matching bonds
            elems = {e0, e1}
            atoms = [(one, two) for one, two in top.bonds
                     if {one.element.symbol, two.element.symbol} == elems]

            # Get indices for the remaining atoms
            indices = []
            for a0, a1 in atoms:
                pair = (a0.index, a1.index)
                # make sure to get the pair in the right order, so that the index
                # for e0 comes before e1
                if a0.element.symbol == e1:
                    pair = pair[::-1]
                indices.append(pair)
            return indices

        nh_donors = get_donors('N', 'H', donor_top)
        oh_donors = get_donors('O', 'H', donor_top)
        xh_donors = np.array(nh_donors + oh_donors)

        if len(xh_donors) == 0:
            # if there are no hydrogens or protein in the trajectory, we get
            # no possible pairs and return nothing
            return np.zeros((0, 3), dtype=int)

        acceptor_elements = frozenset(('O', 'N'))
        acceptors = [a.index for a in acceptor_top.atoms if a.element.symbol in acceptor_elements]

        # Make acceptors a 2-D numpy array
        acceptors = np.array(acceptors)[:, np.newaxis]

        # Generate the cartesian product of the donors and acceptors
        xh_donors_repeated = np.repeat(xh_donors, acceptors.shape[0], axis=0)
        acceptors_tiled = np.tile(acceptors, (xh_donors.shape[0], 1))
        bond_triplets = np.hstack((xh_donors_repeated, acceptors_tiled))

        return bond_triplets

    from mdtraj.geometry.hbond import _compute_bounded_geometry
    distance_cutoff = 0.25  # nanometers
    angle_cutoff = 2.0 * np.pi / 3.0  # radians

    if traj1.topology is None or traj2.topology is None:
        raise ValueError('baker_hubbard requires that traj contain topology '
                         'information')

    # Get the possible donor-hydrogen...acceptor triplets
    h_bond_triplets = _get_bond_triplets_between_structures(traj1.topology, traj2.topology)
    h_bond_triplets_reversed = _get_bond_triplets_between_structures(traj1.topology, traj2.topology)

    mask, distances, angles = _compute_bounded_geometry(traj, bond_triplets,
                                                        distance_cutoff, [1, 2], [0, 1, 2], freq=freq,
                                                        periodic=periodic)

    # Find triplets that meet the criteria
    presence = np.logical_and(distances < distance_cutoff, angles > angle_cutoff)
    mask[mask] = np.mean(presence, axis=0) > freq

    return bond_triplets.compress(mask, axis=0)
