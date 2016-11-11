"""
Truncate: A simple code to truncate DM halos at a given radii.

Input:
------
      String: snaphot name and path e.g: data/halo1.hdf5
      Float: r_cut radii at which to cut the halo.
      String: output file name and path e.g: data/halo1_truncated.hdf5

Usage:
------

python truncate.py snap_name r_cut output
"""

import numpy as np
#from pygadgetreader import readsnap

"""
if (len(sys.argv)!=4):
    print('Error --> Usage: python truncate.py snap_name r_cut output')
    sys.exit()

snap_name = sys.argv[1]
r_cut = float(sys.argv[2])
output = (sys.argv[3])

pos = readsnap(snap_name, 'pos', 'dm')
M = readsnap(snap_name, 'mass', 'dm')
"""

def truncate(pos, vel, rcut, M=False, pot=False):
    truncate = np.where(np.sqrt(pos[:,0]**2.0 + pos[:,1]**2.0 + \
                        pos[:,2]**2.0)<rcut)[0]

    pos_t = pos[truncate]
    vel_t = vel[truncate]
    if  type(M)==type(pos):
        M_t = M[truncate]
        pot_t = pot[truncate]
        return pos_t, vel_t, M_t, pot_t
    else:
        return pos_t, vel_t


"""
f = open(output, "w")
f.write('# posx  posy  posz  Mass')

for i in range(len(truncate)):
    f.write('%f \t %f \t %f \t %f \n'%(posx[i], posy[i], \
            posz[i], M_t[i]))
f.close()
"""
