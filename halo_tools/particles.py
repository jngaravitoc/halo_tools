"""
Code to select host and satellite particles.
The code assumes that the host particles
are first in the ids of the snapshot.

To-Do:
------

Return mass and potential!

"""

import numpy as np



def host_particles(xyz, vxyz, pot=False, mass=False, ids, N_halo):
    """
    Returns the host and the satellite DM particles.

    """
    id_cut = np.sort(ids)[N_halo]
    index_mw = np.where(ids<id_cut)[0]
    index_lmc = np.where(ids>=id_cut)[0]

    return xyz[index_mw], vxyz[index_mw], xyz[index_lmc], vxyz[index_lmc]
