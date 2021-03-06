"""
Script that apply PCA and smoothening
method describes in Weinberg(1996).

Author: Nicolas Garavito-Camargo
08/29/16: University of Arizona.


"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import biff

def read_cov_mat(cov_mat_data, nmax, lmax):
    """
    Function that reads the covariance matrix.
    explain the matrix shape.
    Input:
    ------
    matrix: 1d Array with the covariance matrix
    nmax: int(nmax)
    lmax: int(lmax)

    output:
    -------
    A 2d array with the covariance matrix.

    """
    i=0
    k=0
    j=0

    # Why lmax+2/2?

    covariance_matrix = np.zeros([int((nmax+1)*(lmax+1)*((lmax+2)/2.0)),\
                                  int((nmax+1)*(lmax+1)*((lmax+2)/2.0))])

    for i in range(int((nmax+1)*(lmax+1)*((lmax+2)/2.0))):
        for j in range(int((nmax+1)*(lmax+1)*((lmax+2)/2.0))):
            covariance_matrix[i][j] = cov_mat_data[k]
            k+=1

    return covariance_matrix

def F_k(lambdas):
    F = np.zeros(len(lambdas)+1)
    for i in range(len(lambdas)):
        F[i+1] = sum(lambdas[:i+1])/sum(lambdas)
    return F

def b(a, var_a):

    """
    Function that computing the smoothening
    of the coefficients!

    Input:
    ------

    a: Coefficients.
    var_a: Variance of the coefficients.

    Output:
    -------

    b: The smoothening of the coefficients a.
    """

    b = np.zeros(len(a))
    for i in range(len(b)):
        if var_a[i]==a[i]==0:
            b[i] = 0.0
        else:
            b[i] = 1./(1. + var_a[i]/a[i]**2.0)
    return b

def variance(cov_mat):
    """
    Covariance matrix diagonal.

    Input:
    ------
    cov_matrix: 2d Array with the covariance matrix.

    Output:
    -------
    variance: 1d array with the variance of a.

    """

    var_S = np.zeros(len(cov_mat))
    for i in range(len(cov_mat)):
        var_S[i] = cov_mat[i][i]
    return var_S

def outter_matrix(S):
    """

    """
    out_matrix = np.multiply.outer(S,S)
    lambdas, eigenv = linalg.eig(out_matrix)
    return out_matrix, lambdas, eigenv

def coeff_PCA(S, T):
    new_S = np.dot(T,S).real
    return new_S

def var_PCA(cov_matrix, T_trans):
    T_daga = np.conjugate(T_trans).T
    var_prime = (np.dot(T_trans, np.dot(cov_matrix, T_daga)))
    return np.diagonal(var_prime.real)

def signal2noise(a, var_a):
    SN = (np.abs(a**2.0/var_a))**0.5
    return SN

def mise(b_s, var_a, a):
    D = np.zeros(len(a))
    for i in range(len(D)):
        D[i] = np.nansum(b_s[:i+1]**2.*var_a[:i+1] + \
              (b_s[:i+1]-1)**2.0*a[:i+1]**2.0)
    return D

def nlm_coeff(k, nmax, lmax):

    """
    Returns the coefficients n, l, m for a given number of
    the.....
    e.g k=1 is n=0, l=1, m=0

    Input:
    ------
    k: Number of the coefficient. e.g
    nmax: maximum n in the expansion.
    lmax: maximum l in the expansion.

    Output:
    -------
    (n, l, m)

    """
    i = 0
    for n in range(nmax+1):
        for l in range(lmax+1):
            for m in range(l+1):
                if (i==k):
                    return(n,l,m)
                    break
                i+=1
    print('Error:', k, 'is not in the range of the matrix')


def relevant_coefficients(a, cut, nmax, lmax):

    """
    Find the coefficients that carries more than
    cut% of the information.

    Input:
    ------
    a: Coefficients
    cut: % of
    nmax
    lmax

    Output:
    -------
    3d Matrix ((nmax+1, lmax+1, lmax+1)) with the relevant coefficients.

    """

    coeff_mat = np.zeros((nmax+1, lmax+1, lmax+1))
    all_a = np.abs(a)/np.max(np.abs(a))
    index =  np.where(all_a>cut)[0]
    for i in range(len(index)):
        n, l, m = nlm_coeff(index[i], nmax, lmax)
        coeff_mat[n,l,m] = a[index[i]]
    return coeff_mat

def PCA_coefficients(Snlm, Tnlm, cov_matrix_data, nmax, lmax):

    """
    Compute the SCF coefficients in PCA basis.

    Input:
    ------

    Snlm: 1d Array with the values of Snlm
    Tnlm: 1d Array with the values of Tnlm
    cov_matrix_data: string with the path and name of the covariance
    matrix file.
    nmax: nmax
    lmax: lmax

    Output:
    -------

    Snlm* bnlm*: Snlm coefficient in the principal component basis and
    with smoothening.

    Tnlm* bnlm*: Tnlm coefficient in the principal component basis and
    with smoothening.

    """

    cov_mat_data = np.loadtxt(cov_matrix_data)
    cov_matrix_S = read_cov_mat(cov_mat_data[:,0], nmax, lmax)
    cov_matrix_T = read_cov_mat(cov_mat_data[:,1], nmax, lmax)

    # computing the outter matrix

    outter_mat_S, eigvals_S, Ttrans_S = outter_matrix(Snlm)
    outter_mat_T, eigvals_T, Ttrans_T = outter_matrix(Tnlm)

    # computing the coefficients in the new basis

    Snlm_PCA = coeff_PCA(Snlm, Ttrans_S)
    Tnlm_PCA = coeff_PCA(Tnlm, Ttrans_T)

    # computing the variance in the new basis

    var_S_PCA = var_PCA(cov_matrix_S, Ttrans_S.real)
    var_T_PCA = var_PCA(cov_matrix_T, Ttrans_T.real)

    # computing the smoothening in the PCA basis

    b_S_PCA = b(Snlm_PCA, var_S_PCA)
    b_T_PCA = b(Tnlm_PCA, var_T_PCA)

    new_S = Snlm_PCA*b_S_PCA
    new_T = Tnlm_PCA*b_T_PCA

    return new_S, b_S_PCA, new_T, b_T_PCA


def smooth_coefficients(cov_matrix_data, nmax, lmax, Snlm, Tnlm=False):

    """
    Compute the SCF coefficients in PCA basis.

    Input:
    ------

    Snlm: 1d Array with the values of Snlm
    Tnlm: 1d Array with the values of Tnlm
    cov_matrix_data: string with the path and name of the covariance
    matrix file.
    nmax: nmax
    lmax: lmax

    Output:
    -------

    Snlm* bnlm*: Snlm smoothed coefficients.

    Tnlm* bnlm*: Tnlm smoothed coefficients.

    """
    cov_mat_data = np.loadtxt(cov_matrix_data)

    if Tnlm[0]!=-99999:
        cov_matrix_T = read_cov_mat(cov_mat_data[:,1], nmax, lmax)
        var_T = np.diagonal(cov_matrix_T)
        b_T = b(Tnlm, (var_T))
        new_T = Tnlm*b_T


        cov_mat_data = np.loadtxt(cov_matrix_data)
        cov_matrix_S = read_cov_mat(cov_mat_data[:,0], nmax, lmax)
        var_S = np.diagonal(cov_matrix_S)
        b_S = b(Snlm, (var_S))
        new_S = Snlm*b_S #/np.max(np.abs(b_S))
        return new_S, new_T

    else:
        cov_matrix_S = read_cov_mat(cov_mat_data, nmax, lmax)
        var_S = np.diagonal(cov_matrix_S)
        b_S = b(Snlm, (var_S))
        new_S = Snlm*b_S #/np.max(np.abs(b_S))

    return new_S

def pot_scatter(potential, posx, posy, title, xmin=False, xmax=False, \
                ymin=False, ymax=False, sf=False, figt=False):
    """
    Plot potential contours of a distribution of particles with given
    potential.

    Input:
    ------
    potential: 1d array of the potential of the particles in the halo.

    posx: 1d array with the positions of the particles in an axys to
    project the potential.

    posy: 1d array with the positions of the particles in an axys to
    project the potential.

    title: plot title

    ouput:
    ------
    contour plot.

    """

    # Picking colors
    N_cuts = 20
    colors = np.r_[np.linspace(0.1, 1, N_cuts), np.linspace(0.1, 1, N_cuts)]
    cm = plt.get_cmap('Spectral')
    my_colors = cm(colors)

    # Potential levels
    pot_cuts_all_nb = np.linspace(min(np.abs(potential)), \
                      max(np.abs(potential)), N_cuts)

    plt.figure(figsize=(5,5))
    plt.text(-500, 500, ('%.2e'%(min(np.abs(potential)))))
    plt.text(-500, 420, ('%.2e'%(max(np.abs(potential)))))

    plt.title(title, fontsize=30)
    for i in range(1,N_cuts):
        index_c = np.where((np.abs(potential)<pot_cuts_all_nb[i]) &
                           (np.abs(potential)>pot_cuts_all_nb[i-1]))[0]

        plt.scatter(posx[index_c], posy[index_c], c=my_colors[i],\
                    edgecolors='none', s=1)


    plt.xlabel('$Y[kpc]$', fontsize=25)
    plt.ylabel('$Z[kpc]$', fontsize=25)
    if xmin:
        plt.xlim(xmin, xmax)
        plt.ylim(ymax, ymin)
    if sf:
        plt.savefig(figt+'.png', bbox_inches='tight', dpi=300)

# plot potenital contours with no scatter plot.


def plot_potential_contours(idx, S, T, grid_max=250., n_grid=300):
    grid = np.zeros((3,n_grid*n_grid))
    grid[idx] = np.vstack(map(np.ravel, np.meshgrid(\
                           np.linspace(-grid_max,grid_max,n_grid),\
                           np.linspace(-grid_max,grid_max,n_grid))))

    pot = biff.potential(np.ascontiguousarray(grid.T), S, T, M=1.0, r_s=40.85)

    pot=np.abs(pot) 
    #think in how to plot negative potentials in the logscale

    percs = np.log10([pot[pot > 0].min(), pot[pot > 0].max()])
    #percs = ([dens.max(), dens.min()])

    levels = np.logspace(percs[0], percs[1], 20)

    shp = (n_grid, n_grid)

    fig,ax = plt.subplots(1,1,figsize=(8,6))


    cf0 = ax.contourf(grid[idx[0]].reshape(shp), grid[idx[1]].reshape(shp),\
                      pot.reshape(shp), cmap='Spectral', levels=levels,extend='both')

    cbar0 = plt.colorbar(cf0,)

    ax.set_xlim(-grid_max, grid_max)
    ax.set_ylim(-grid_max, grid_max)

    ax.set_xlabel("${}$".format('xyz'[idx[0]]), fontsize=25)
    ax.set_ylabel("${}$".format('xyz'[idx[1]]), fontsize=25)

    tmp = np.array([0,1,2])
    tmp = np.delete(tmp, np.where((tmp==idx[0]) | (tmp == idx[1]))[0])
    #ax.set_title("Isodensity contours at${}=0$".format('xyz'[tmp[0]]))

    fig.tight_layout()

    return fig,ax

# function to plot multipole moments pontential




