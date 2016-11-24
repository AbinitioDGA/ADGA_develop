#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
#from mpl_toolkit import Axes3D

# create
nx=20
ny=14
nz=10
nk=nx*ny*nz
kmeshx = np.linspace(0, 1, nx, endpoint=False)
kmeshy = np.linspace(0, 1, ny, endpoint=False)
kmeshz = np.linspace(0, 1, nz, endpoint=False)
t=0.5/np.sqrt(6.0)
Hk = -2.0*t*(np.cos(kmeshz*2*np.pi)[None,None,:] + np.cos(kmeshy*2*np.pi)[None,:,None] + np.cos(kmeshx*2*np.pi)[:,None,None])

# write
# k nat nd np nlig
print "%d 1 1 0 0" % nk
for pos, H in np.ndenumerate(Hk):
    print "%g\t%g\t%g" % (kmeshx[pos[0]], kmeshy[pos[1]], kmeshz[pos[2]])
    print " %g\t0.0" % H

