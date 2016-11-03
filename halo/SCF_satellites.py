import numpy as np
from pygadgetreader import readsnap
from .shapes import *
from .particles import *
import biff



def recenter_host(pos, centers):
    """
    Recenter a halo in its cm.

    Parameters:
    -----------
    pos: ndarray [Npart, pos]
         Array with the positions of the halo.
    centers: array [cmx, cmy, cmz]
         Array with the position of the halo centers

    """
    pos[0] = pos[0] - centers[0]
    pos[1] = pos[1] - centers[1]
    pos[2] = pos[2] - centers[2]

    return pos

def reading_snaps(snap_path, snap_name, N_snap, N_part_host):
    """
    function that read a snapshot using pygadget reader
    """
    pos = readsnap(lmc_path + lmc_snap + '_{:03d}'.format(N_snap), 'pos', 'dm')
    M = readsnap(lmc_path + lmc_snap + '_{:03d}'.format(N_snap),'mass', 'dm')
    vel = readsnap(lmc_path + lmc_snap + '_{:03d}'.format(N_snap),'vel', 'dm')
    pot = readsnap(lmc_path + lmc_snap + '_{:03d}'.format(N_snap),'pot', 'dm')
    pids = readsnap(lmc_path + lmc_snap + '_{:03d}'.format(N_snap),'pid', 'dm')

    MW_pos, MW_vel, MW_mass, MW_pot, LMC_pos, LMC_vel, LMC_mass, \
    LMC_pot = host_sat_prticles(pos, vel, pids,  N_part_host, pot, M)

    return MW_pos, MW_vel, MW_mass, MW_pot, LMC_pos, LMC_vel, \
           LMC_mass, LMC_pot

def trunc_radii_sat(pos, rcut):

    """
    Find where the halo stops being spherical. By computing
    the eigenvalues of the 'Inertia Tensor'.

    Parameters:
    -----------
    pos: ndarray
         Halo positions of the halo

    rcut: float
         Maximum radii to look for the spherical radius.

    Returns:
    --------
    rcut: float
          max radius at which the halo is still spherical.
    """

    r = np.linspace(1, rcut, 100)
    s=1
    q=1
    i=0
    while ((np.abs(s-1)< 0.1) & (np.abs(q-1)<0.1) & (i<200)):
        s, q = shapes.iterate_volume(pos[:,0], pos[:,1], pos[:,2], r[i], 1E-1)
        i+=1
    return r[i-2]


def mask(pos, M, pot, center, r_cut):
    """

    """
    inn = np.where(np.sqrt((pos[:,0]-xcm)**2.0 + (pos[:,1]-ycm)**2.0 + (pos[:,2]-zcm)**2.0)<r_cut)[0]
    out = np.where(np.sqrt((pos[:,0]-xcm)**2.0 +(pos[:,1]-ycm)**2.0 + (pos[:,2]-zcm)**2.0)>=r_cut)[0]
    return pos[inn], M[inn], pot[inn], pos[out], M[out], pot[out]


def mask_satellite(snap_path, snap_name, orbits, Ni, Nf, rcut):
    for i in range(Ni, Nf+1):
        MW_pos, MW_vel, MW_mass, MW_pot, LMC_pos, LMC_vel, LMC_mass, \
        LMC_pot = reading_snaps(snap_path, snap_name, N_snap, N_part_host)

        sat_centerx = orbits[i,5]
        sat_centery = orbits[i,6]
        sat_centerz = orbits[i,7]

        LMC_pos_c = recenter_host(LMC_pos, [sat_centerx, sat_centery, sat_centerz])

        r_cut = trunc_radii_sat(LMC_pos_c, rcut)

        sat_pos_inn, sat_M_inn, sat_pot_inn, sat_pos_out, sat_M_out,\
        sat_pot_out = mask(LMC_pos, lmc_orbit[i,2], lmc_orbit[i,3],\
                           lmc_orbit[i,4], r_cut, LMC_mass, LMC_pot)

        #host_pos_inn, host_M_inn, host_pot_inn, host_pos_out, \
        #host_M_out, host_pot_out = mask(MW_pos, lmc_orbit[i,2],\
        #                                lmc_orbit[i,3], lmc_orbit[i,4],\
        #                                r_cut, MW_mass, MW_pot)

        out_pos = np.concatenate((MW_pos, sat_pos_out), axis=0)
        out_M = np.concatenate((MW_mass, sat_M_out), axis=0)
        out_pot = np.concatenate((MW_pot, sat_pot_out), axis=0)

        return out_pos, out_M, out_pot, sat_pos_inn, sat_M_inn\
               ,sat_pot_inn

