import numpy as np
import biff


def regions(mw_inn, mw_out, lmc_inn, lmc_out):
    lmc_posx = lmc_inn[:,0]
    lmc_posy = lmc_inn[:,1]
    lmc_posz = lmc_inn[:,2]
    lmc_pot = lmc_inn[:,3]
    lmc_m = lmc_inn[:,4]

    return lmc_posx, lmc_posy
