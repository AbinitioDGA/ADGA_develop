#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import h5py
import optparse
import sys

def componentband2index(i1,i2,i3,i4,ndim):
    return ndim**3*(i1-1) + ndim**2*(i2-1) + ndim*(i3-1) + i4

parser = optparse.OptionParser()
parser.add_option("--ndim", type="int", dest="ndim")
options, args = parser.parse_args()

if len(args) != 2:
    parser.error("expected input text file and output hdf5 file")

infile  = args[0]
outfile = args[1]

ndim    = options.ndim

# open the input text file and count the lines
f = open(str(infile), 'r')

try:
    g = h5py.File(str(outfile), 'w-')
except IOError:
    print('Output hdf5 file already exists, cannot overwrite.')
    print('Remove existing file or choose another name.')
    sys.exit()

# read R-points and create hdf5 group
Rpoints = np.genfromtxt(str(infile), dtype=np.float64, usecols=(0,1,2)) # always 3-dimensional
g.create_group('.axes')
g.create_dataset('.axes/R-points', data=Rpoints, dtype=np.float64)

# the output format from Jan's U-files
# 1 1 x x ... x = 1, ndim
# (cont) 1 x x 1 ... x = 1, ndim
# (cont) 1 x 1 x ... x = 1, ndim

# Jan's format conforms the following diagram: a-b-c-d
#  a    _________________  d
#           |       |
#           |       |
#           |       |
#  b    ____|_______|____  c
#
#
# our format (ADGA) conforms the following diagram: a-b-c-d
#  a    _________________  d
#           |       |
#           |       |
#           |       |
#  c    ____|_______|____  b
#
# this means we have to switch the second and third leg

for i in xrange(1,ndim+1):
    # we switch from 1 1 x x to 1 x 1 x
    ind = componentband2index(1,i,1,i,ndim)
    column = 2 + i # first three columns are the r-coordinates
    Rdata = np.genfromtxt(str(infile), dtype=np.float64, usecols=(column, ))
    print(Rdata[:5])
    g.create_dataset('{:05}/value'.format(ind), data=Rdata)

for i in xrange(2,ndim+1): # we ignore the first column of this subsection because it is 1111 which we already got from the part above
    # we switch from 1 x x 1 to 1 x x 1
    ind = componentband2index(1,i,i,1,ndim)
    column = 2 + ndim + i
    Rdata = np.genfromtxt(str(infile), dtype=np.float64, usecols=(column, ))
    print(Rdata[:5])
    g.create_dataset('{:05}/value'.format(ind), data=Rdata)

for i in xrange(2,ndim+1): # same here
    # we switch from 1 x 1 x to 1 1 x x
    ind = componentband2index(1,1,i,i,ndim)
    column = 2 + ndim + ndim + i
    Rdata = np.genfromtxt(str(infile), dtype=np.float64, usecols=(column, ))
    print(Rdata[:5])
    g.create_dataset('{:05}/value'.format(ind), data=Rdata)

sys.exit()

# now we use the orbital symmetry

if ndim > 1:
    ind_from = componentband2index(1,1,1,1,ndim)
    for i in xrange(2,ndim+1):
        ind_to   = componentband2index(i,i,i,i,ndim)
        g['{:05}/value'.format(ind_to)] = g['{:05}/value'.format(ind_from)] # h5py HardLink

    ind_from = componentband2index(1,2,1,2,ndim)
    for i in xrange(1,ndim+1):
        for j in xrange(1,ndim+1):
            if i == j:
                continue
            ind_to   = componentband2index(j,i,j,i,ndim)
            if not '/{:05}'.format(ind_to) in g:
                g['{:05}/value'.format(ind_to)] = g['{:05}/value'.format(ind_from)] # h5py HardLink

    ind_from = componentband2index(1,2,2,1,ndim)
    for i in xrange(1,ndim+1):
        for j in xrange(1,ndim+1):
            if i == j:
                continue
            ind_to   = componentband2index(j,i,i,j,ndim)
            if not '/{:05}'.format(ind_to) in g:
                g['{:05}/value'.format(ind_to)] = g['{:05}/value'.format(ind_from)] # h5py HardLink

    ind_from = componentband2index(1,1,2,2,ndim)
    for i in xrange(1,ndim+1):
        for j in xrange(1,ndim+1):
            if i == j:
                continue
            ind_to   = componentband2index(i,i,j,j,ndim)
            if not '/{:05}'.format(ind_to) in g:
                g['{:05}/value'.format(ind_to)] = g['{:05}/value'.format(ind_from)] # h5py HardLink


g.close()
print('Done.')
