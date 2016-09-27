#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Code to compute the center of Mass of DM halos.

Author: J. Nicolas Garavito-Camargo.

09/22/16

University of Arizona.

To-Do:
------
1. Work in the main, writing orbit function.
2. Verify the Shrinking Sphere algorithm.

"""

import numpy as np
from pygadgetreader import *
import sys
import argparse

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

    return np.array([xcm, ycm, zcm]), np.array([vxcm, vycm, vzcm])


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


def CM_iterative(xyz, vxyz, delta=0.025):

    """
    Compute the center of mass coordinates and velocities of a halo
    using the Shrinking Sphere Method Power et al 2003.
    It iterates in radii until reach a convergence given by delta
    or 1% of the total number of particles.

    Parameters:
    -----------
    xyz: numpy.ndarray (n,3)
      cartesian coordinates of the halo

    vxys: numpy.ndarray (n,3)
      cartesian velocities of the halo

    delta(optional): float (defautl=0.025kpc)
      Precision of the CM,

    Returns:
    --------
    rcm: numpy.ndarray (1,3
      cartesian coordinates of the CM.
    vcm: numpy.ndarray(1,3)
      cartesian velocity the CM.

    """
    rCM_new, vCM_new = CM_basic(xyz, vxyz)

    xCM = 0.0
    yCM = 0.0
    zCM = 0.0


    while (np.sqrt((rCM_new[0]-xCM)**2 + (rCM_new[1]-yCM)**2 + (rCM_new[2]-zCM)**2) > delta):

        xCM = rCM_new[0]
        yCM = rCM_new[1]
        zCM = rCM_new[2]

        R = np.sqrt((xyz[:,0] - xCM)**2 + (xyz[:,1] - yCM)**2 + (xyz[:,2] - zCM)**2)
        Rmax = max(R)

        # Rate at which the sphere is reducing with the enclosed
        # particles.reducing 2.5%
        index = np.where(R<Rmax*0.975)[0]
        xyz = xyz[index]
        vxyz = vxyz[index]

        rCM_new, vCM_new = CM_basic(xyz, vxyz)

    return rCM_new, vCM_new

# function that computes the CM using the 10% most bound particles!
# I am using the potential method any more, its not useful to find the 
#LMC CM because the particles feel the potential of the MW.
"""
def potential_CM(potential, x, y, z, vx, vy, vz):
    index = np.where(potential< min(potential)*0.90)[0]
    x_p = x[index]
    y_p = y[index]
    z_p = z[index]
    vx_p = vx[index]
    vy_p = vy[index]
    vz_p = vz[index]
    N = len(x_p)
    x_cm = sum(x_p)/N
    y_cm = sum(y_p)/N
    z_cm = sum(z_p)/N
    vx_cm = sum(vx_p)/N
    vy_cm = sum(vy_p)/N
    vz_cm = sum(vz_p)/N
    Rcm = np.sqrt(x_cm**2.0 + y_cm**2.0 + z_cm**2.0)
    Vcm = np.sqrt(vx_cm**2.0 + vy_cm**2.0 + vz_cm**2.0)
    return x_cm, y_cm, z_cm, vx_cm, vy_cm, vz_cm, Rcm, Vcm


#Function that computes the CM of the halo using the minimum of the
#potential:

def CM_host_potential(xyz, vxyz, pot, Rinn=2.0):

"""

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
            MW_xyz, MW_vxyz, LMC_xyz, LMC_vxyz = MW_LMC_particles(xyz, vxyz, pids, NMW_particles)
            MW_rcm[i], MW_vcm[i] = CM(MW_xyz, MW_vxyz, delta)
            LMC_rcm[i], LMC_vcm[i] = CM(LMC_xyz, LMC_vxyz, delta)
        else:
            MW_rcm[i], MW_vcm[i] = CM(xyz, vxyz, delta)
    return MW_rcm, MW_vcm, LMC_rcm, LMC_vcm


def writing_file(out_name, Host_rcm, Host_vcm, Sat_rcm ,Sat_vcm):
    """

    Function that writes the output.

    """
    f = open(out_name, 'w')
    f.write("#Time(Gyrs) | Rgal(kpc) | Xsat[kpc] | Ysat[kpc] | Zsat[kpc] |Xhost[kpc] | Yhost[kpc] Zhost[kpc] |"\
            "Vgal | Vxsat | Vysat | Vzsat | Vxhost | Vyhost | Vzhost |\n")

    for i in range(0, len(Rgal)):
        f.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%(Rgal[i], Xsat[i], Ysat[i],\
        Zsat[i], X[i], Y[i], Z[i], Vgal[i], VXsat[i], VYsat[i], VZsat[i], VX[i], VY[i], VZ[i]))
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
