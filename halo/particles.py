"""
Code to select host and satellite particles.
The code assumes that the host particles
are first in the ids of the snapshot.

To-Do:
------

Return mass and potential!

"""

import numpy as np

def host_sat_particles(xyz, vxyz, ids, N_halo, pot=False, mass=False):
    """
    Returns the host and the satellite DM particles.

    """
    id_cut = np.sort(ids)[N_halo]
    index_mw = np.where(ids<id_cut)[0]
    index_lmc = np.where(ids>=id_cut)[0]

    # Look for the proper way to do this!
    if type(pot)==type(xyz):
        return xyz[index_mw], vxyz[index_mw], mass[index_mw],\
               pot[index_mw], xyz[index_lmc], vxyz[index_lmc],\
               mass[index_lmc], pot[index_lmc]
    else:
        return xyz[index_mw], vxyz[index_mw], xyz[index_lmc], vxyz[index_lmc]
