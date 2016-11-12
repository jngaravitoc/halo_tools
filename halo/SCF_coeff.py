import numpy as np
import biff
from pygadgetreader import readsnap
from scipy import interpolate

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

def read_coeff(path, snap_name, i):
    coeffs = np.loadtxt(path+snap_name+'{:d}'.format(i)+'_n10l5.txt',skiprows=1)
    S_host = coeffs[:,0]
    T_host = coeffs[:,1]
    S_sat = coeffs[:,2]
    T_sat = coeffs[:,3]
    return S_host, T_host, S_sat, T_sat

def ST_2interpol(path, snap_name, nmax, lmax, n_snap, h, dt_snap):
    #n_snap = 8
    #nmax = 1
    #lmax = 1
    #dt_snap = 0.02
    #h = 0.001

    # Number of coefficients Snlm / Tnlm
    n_coeff = ((nmax+1)*(lmax+1)*(lmax+1))  # ** this might change orbe general

    # S, T matrices, rows: time / cols: Number of coefficients
    S_host_matrix = np.zeros((n_snap, n_coeff))
    T_host_matrix = np.zeros((n_snap, n_coeff))
    S_sat_matrix = np.zeros((n_snap, n_coeff))
    T_sat_matrix = np.zeros((n_snap, n_coeff))

    # Finish this orbit upload!!
    #rcm_mw, vcm_mw, rcm_lmc, vcm_lmc = np.loadtxt('../data/orbits/')
    # Filling the matrix with the coefficients form the N-body sim.
    for i in range(n_snap):
        S_host_matrix[i], T_host_matrix[i], S_sat_matrix[i], T_sat_matrix[i] = read_coeff(path, snap_name, i)

    # time: Time between every snapshot
    time = np.linspace(0, dt_snap*(n_snap-1), n_snap)

    # N_interp: Number of times in the interpolation (this depent in
    # the resolution of the leapfrog integration h)
    N_interp = (dt_snap*(n_snap-1)/h)+1

    # t: Time between every time step in the interpolation
    t = np.linspace(0, dt_snap*(n_snap-1), N_interp)

    #S_matrix_interp: Matrix with the new interpolated Snlm (times, Svalues)
    S_host_matrix_interp = np.zeros((len(t), n_coeff))
    T_host_matrix_interp = np.zeros((len(t), n_coeff))
    S_sat_matrix_interp = np.zeros((len(t), n_coeff))
    T_sat_matrix_interp = np.zeros((len(t), n_coeff))

    # Interpolating 
    for i in range(n_coeff):
        fs_host = interpolate.interp1d(time, S_host_matrix[:,i], kind='linear')
        ft_host = interpolate.interp1d(time, T_host_matrix[:,i], kind='linear')
        fs_sat = interpolate.interp1d(time, S_sat_matrix[:,i], kind='linear')
        ft_sat = interpolate.interp1d(time, T_sat_matrix[:,i], kind='linear')
        S_host_matrix_interp[:,i]= fs_host(t)
        T_host_matrix_interp[:,i]= ft_host(t)
        S_sat_matrix_interp[:,i]= fs_sat(t)
        T_sat_matrix_interp[:,i]= ft_sat(t)

    S_host_new = np.reshape(S_host_matrix_interp, (N_interp, nmax+1, lmax+1, lmax+1))
    T_host_new = np.reshape(T_host_matrix_interp, (N_interp, nmax+1, lmax+1, lmax+1))
    S_sat_new = np.reshape(S_sat_matrix_interp, (N_interp, nmax+1, lmax+1, lmax+1))
    T_sat_new = np.reshape(T_sat_matrix_interp, (N_interp, nmax+1, lmax+1, lmax+1))

    return S_host_new, T_host_new, S_sat_new, T_sat_new
