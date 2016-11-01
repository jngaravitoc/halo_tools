#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Code to compute the center of Mass of DM halos.

Author: J. Nicolas Garavito-Camargo.

09/22/16

University of Arizona.

To-Do:
------
1. Implement the main function.

"""

import numpy as np
from pygadgetreader import *
from .particles import host_sat_particles
from .truncate import truncate

def CM_basic(xyz, vxyz):

    """
    Computing the CM of particles using the standard equation:
    r_cm = sum(r_i) / N

    Parameters:
    -----------
    xyz: numpy.ndarray
      Positions of the particles

    vxyz: numpy.ndarray
      Velocities of the particles

    Output:
    -------
    r_cm: numpy.ndarray
      coordinates of the CM (xcm, ycm, zcm)
    vr_cm numpy.ndarray
      velocity of the CM (vxcm, vycm, vzcm)

    """
    assert (len(xyz[:,0])==len(xyz[:,1]) == len(xyz[:,2])== \
            len(vxyz[:,0])==len(vxyz[:,1])==len(vxyz[:,2]))

    N = len(xyz[:,0])

    xcm = sum(xyz[:,0])/N
    ycm = sum(xyz[:,1])/N
    zcm = sum(xyz[:,2])/N

    vxcm = sum(vxyz[:,0])/N
    vycm = sum(vxyz[:,1])/N
    vzcm = sum(vxyz[:,2])/N

    return np.array([xcm, ycm, zcm]).T, np.array([vxcm, vycm, vzcm]).T


def CM_disk_potential(xyz, vxyz, pot, Inn_radius=2.0):

    """
    CM_disk_potential:
    ------------------

    Computes the CM of disk particles by finding the CM of the
    most bound particles. i.e particles with the minimum value
    of the gravitational potential.

    Parameters:
    -----------

    xyz: numpy.ndarray
      Particle positions.
    vxyz: numpy.ndarray
      Particle velocities.
    pot: numpy.ndarray
      Particle potential.
    Inn_radius: float, default 2kpc
      Inner radius at which the most bound
      particles are going to be.

    Output:
    -------
    r_cm: numpy.ndarray
      positions of the CM
    v_cm: numpy.ndarray
      velocity of the CM

    """

    min_pot = np.where(pot==min(pot))[0]
    x_min = xyz[min_pot,0]
    y_min = xyz[min_pot,1]
    z_min = xyz[min_pot,2]

    avg_particles = np.where(np.sqrt((xyz[:,0]-x_min)**2.0 + (xyz[:,1]-y_min)**2.0 +\
                            (xyz[:,2]-z_min)**2.0)<Inn_radius)[0]

    r_cm, v_cm = CM_basic(xyz[avg_particles], vxyz[avg_particles])

    return r_cm, v_cm


def shrinking_sphere(xyz, vxyz, delta=0.025):

    """

    Compute the center of mass coordinates and velocities of a halo
    using the Shrinking Sphere Method Power et al 2003.
    It iterates in radii until reach a convergence given by delta
    or 1% of the total number of particles.

    Parameters:
    -----------
    xyz: numpy.ndarray (n,3)
      Cartesian coordinates of the halo

    vxys: numpy.ndarray (n,3)
      Cartesian velocities of the halo

    delta(optional): float (default=0.025kpc)
      Precision of the CM,

    Returns:
    --------
    rcm: numpy.ndarray (1,3
      Cartesian coordinates of the CM.
    vcm: numpy.ndarray(1,3)
      Cartesian velocity the CM.

    """
    rCM_new, vCM_new = CM_basic(xyz, vxyz)

    xCM = 0.0
    yCM = 0.0
    zCM = 0.0

    Ni = len(xyz)

    while ((np.sqrt((rCM_new[0]-xCM)**2 + (rCM_new[1]-yCM)**2 + (rCM_new[2]-zCM)**2) > delta)
           & (len(xyz)>0.01*Ni)):

        xCM = rCM_new[0]
        yCM = rCM_new[1]
        zCM = rCM_new[2]

        R = np.sqrt((xyz[:,0] - xCM)**2 + (xyz[:,1] - yCM)**2 + (xyz[:,2] - zCM)**2)
        Rmax = max(R)

        # Rate at which the sphere is reducing with the enclosed
        # particles reducing radii = 2.5%
        index = np.where(R<Rmax*0.5)[0]
        xyz = xyz[index]
        vxyz = vxyz[index]

        rCM_new, vCM_new = CM_basic(xyz, vxyz)

    return rCM_new, vCM_new



#Function that computes the CM of the halo using the minimum of the
#potential:



def orbit_CM(path, snap_name, initial_snap, final_snap, NMW_particles,\
             delta=0.025, lmc=False):

    """
    Computes the orbit of the host and the satellite.
    The COM of the host and the satellite are compute using the
    shrinking sphere method at each snapshot.

    Parameters:
    -----------

    path: Path to the simulation snapshots
    snap_name: Base name of the snaphot without the number and
    file type, e.g: LMCMW
    initial_snap: Number of the initial snapshot
    final_snap: Number of the final snapshot
    NMW_particles: Number of MW particles in the simulation.
    delta: convergence distance
    lmc: track the lmc orbit. (default = False)

    Output:
    -------
    XMWcm, vMWcm, xLMCcm, vLMCcm: 4 arrays containing the coordinates
    and velocities of the center of mass with reference to a (0,0,0) point
    at a given time.

    """
    r_mwcut = 500
    r_lmccut = 100
    N_snaps = final_snap - initial_snap + 1
    MW_rcm = np.zeros((N_snaps,3))
    MW_vcm = np.zeros((N_snaps,3))
    LMC_rcm = np.zeros((N_snaps,3))
    LMC_vcm = np.zeros((N_snaps,3))

    for i in range(initial_snap, final_snap+1):
        xyz = readsnap(path + snap_name + '_{:03d}.hdf5'.format(i),'pos', 'dm')
        vxyz = readsnap(path + snap_name +'_{:03d}.hdf5'.format(i),'vel', 'dm')
        pids = readsnap(path + snap_name +'_{:03d}.hdf5'.format(i),'pid', 'dm')
        if lmc==True:
            MW_xyz, MW_vxyz, LMC_xyz, LMC_vxyz = host_sat_particles(xyz, vxyz, pids, NMW_particles)
            MW_pos_t, MW_vel_t = truncate(MW_xyz, MW_vxyz, r_mwcut)
            LMC_pos_t, LMC_vel_t = truncate(LMC_xyz, LMC_vxyz, r_lmccut)
            MW_rcm[i-initial_snap], MW_vcm[i-initial_snap] = shrinking_sphere(MW_pos_t, MW_vel_t, delta)
            LMC_rcm[i-initial_snap], LMC_vcm[i-initial_snap] = shrinking_sphere(LMC_pos_t, LMC_vel_t, delta)
        else:
            MW_pos_t, MW_vel_t = truncate(xyz, vxyz, r_mwcut)
            MW_rcm[i-initial_snap], MW_vcm[i-initial_snap] = shrinking_sphere(MW_pos_t, MW_vel_t, delta)

    return MW_rcm, MW_vcm, LMC_rcm, LMC_vcm


def writing_file(out_name, time, Host_rcm, Host_vcm, Sat_rcm ,Sat_vcm):
    """

    Function that writes the output.

    """

    MW_r = np.sqrt(Host_rcm[0]**2.0 + Host_rcm[1]**2.0 + Host_rcm[2]**2.0)
    MW_v = np.sqrt(Host_vcm[0]**2.0 + Host_vcm[1]**2.0 + Host_vcm[2]**2.0)
    R_gal = np.sqrt((Host_rcm[0]-Sat_rcm[0])**2.0 + (Host_rcm[1]-\
                    Sat_rcm[1])**2.0 +(Host_rcm[2]-Sat_rcm[2])**2.0)


    V_gal = np.sqrt((Host_vcm[0]-Sat_vcm[0])**2.0 + (Host_vcm[1]-\
                    Sat_vcm[1])**2.0 +(Host_vcm[2]-Sat_vcm[2])**2.0)
    f = open(out_name, 'w')
    f.write("#Time (Gyrs) | R_mw(kpc) | Rgal_sat(kpc) |  V_mw(km/s) "\
            "| V_sat(km/s) |  Xhost[kpc] | Yhost[kpc] | Zhost[kpc] |"\
            " Xsat[kpc] | Ysat[kpc] | Zsat[kpc] | Vxhost | Vyhost |"\
            "Vzhost | Vxsat | Vysat | Vzsat |\n")

    for i in range(0, len(R_gal)):
        f.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"\
                %(time[i], MW_r[i], R_gal[i], MW_v[i], V_gal[i],\
                Host_rcm[i,0], Host_rcm[i,1], Host_rcm[i,2], Sat_rcm[i,0], Sat_rcm[i,1],\
                Sat_rcm[i,2], Host_vcml[i,0], Host_vcm[i,1],Host_vcm[i,2], Sat_vcm[i,0],\
                Sat_vcm[i,1], Sat_vcm[i,2]))
    f.close()

if __name__ == '__main__':

    if len(sys.argv) != 8:
        print('Usage: python orbit_cm.py snap_base_name inital_snap_number'\
              'final_snap_number path2file out?name  #DMhost #DMsat')
        print('Ex: python orbit_cm.py snap 0 50 pat2file out_name Nhost Nsat')
        sys.exit(0)


    # global variables to read the snapshots
    snap = str(sys.argv[1])
    # Initial and final snapshot number
    i_n = int(sys.argv[2])
    i_f = int(sys.argv[3])
    #Output Name
    out_name = str(sys.argv[5])
    #Number of particles of the host
    Nhost = int(sys.argv[6])
    #Number of particles of the sat
    Nsat = int(sys.argv[7])
    path = str(sys.argv[4]) #'../../data/LMCMW/MW1LMC4/a1/'

    # Number of Snapshots
    MW_rcm, MW_vcm, LMC_rcm, LMC_vcm = orbit_CM(path, snap, i_n, i_f, Nhost)

    D = 0.025 # precision of the CM computation in Kpc

    print('Writing the data')
