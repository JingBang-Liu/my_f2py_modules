import numpy as np
import os
import readBins


path = os.getcwd()
nbin = np.array((153,4,8))
timemarks = 500
filename = path + "/../../MD/MD_results/rupture_1/1/bins_dens.dat"

a = readBins.readbins(timemarks,153,4,8,filename)
a = np.array(a)
print(a.shape)
