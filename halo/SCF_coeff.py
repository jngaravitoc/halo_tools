import numpy as np
import biff
from pygadgetreader import readsnap
#import particles
#import truncate

def clean_halo(path, snap_name, N):
    pos = readsnap(path + snap_name + '_{:03d}'.format(N), 'pos', 'dm')
    M = readsnap(path + snap_name + '_{:03d}'.format(N), 'mass', 'dm')
    vel = readsnap(path + snap_name + '_{:03d}'.format(N), 'vel','dm')
    pot = readsnap(path + snap_name + '_{:03d}'.format(N), 'pot','dm')
    pids = readsnap(path + snap_name + '_{:03d}'.format(N), 'pid','dm')

    particles
    cm
    truncate
    return pos,  M

def compute_coeff(path, snap_name, N, lmax, nmax):
    pos, mass = clean_halo(path, snap_name, N)
    S, T = biff.compute_coeffs_discrete(pos.astype(np.float64), mass.astype(np.float64), nmax, lmax, r_s)
    return S, T

def write_coeff(path, snap_name, N, lmax, nmax, file_name):
    S, T = compute_coeff(path, snap_name, N, lmax, nmax)
    f = open(file_name, 'w')
    S_flat = S.flatten()
    T_flat = T.flatten()
    for i in range(len(S_flat)):
        f.write("%f \t %f \n"%(S_flat[i], T_flat[i]))
    f.close()


#def read_coeff(nmax, lmax):
#    f = read
#    S = reshape().((nmax, lmax, lmax))

