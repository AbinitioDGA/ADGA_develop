#!/usr/bin/env python

# This script produces an Abinitio DGA input file from
# a two-particle Green's function that was measured 
# with segment code Z sampling.

# It is necessary to do the following transposition:
# 
# G(adga)v,v',w,_abcd = G(seg)v',v,-w,_cdab
# 
# The transposition in the spin/band indices 
# can be omitted in SU(2) symmetric 1-band case.

import numpy as np
import h5py
import optparse


parser = optparse.OptionParser(usage = "%prog [OPTIONS] <hdf file> ...")


(options, args) = parser.parse_args()

if not args:
  parser.error("expecting two HDF5 filenames and beta")

segfilename=args[0]
adgafilename=args[1]


segfile=h5py.File(segfilename,'r')
adgafile=h5py.File(adgafilename,'w')

nineq = segfile['.config'].attrs['general.nat']
beta  = segfile['.config'].attrs['general.beta']
boslist  = segfile['.axes/iwb-g4'].value # 2.*np.pi/beta*np.arange(-nb,nb+1)
fermlist = segfile['.axes/iwf-g4'].value # np.pi/beta*(2.*np.arange(-nf,nf)+1)

adgafile.create_group('.axes')
adgafile['.axes/iwb-g4']=boslist
adgafile['.axes/iwf-g4']=fermlist

for ineq in xrange(nineq):
  g4seg=segfile['stat-001/ineq-{:03}/g4iw/value'.format(ineq+1)].value

  g4dga=g4seg.transpose((0,1,2,3,5,4,6))[...,::-1]

  nb=g4dga.shape[-1]//2
  nf=g4dga.shape[-2]//2
  adgafile.create_group('ineq-{:03}'.format(ineq+1))
  adgafile.create_group('ineq-{:03}/dens'.format(ineq+1))
  adgafile.create_group('ineq-{:03}/magn'.format(ineq+1))

  # channel symmetrization
  g4dens=0.5*(g4dga[0,0,0,0]+g4dga[0,1,0,1]+g4dga[0,0,0,1]+g4dga[0,1,0,0])
  g4magn=0.5*(g4dga[0,0,0,0]+g4dga[0,1,0,1]-g4dga[0,0,0,1]-g4dga[0,1,0,0])


  # create groups for the actual data.
  # structure: ineq/chann/bosfreq/bandgroup
  for i in xrange(2*nb+1):
    adgafile['ineq-{:03}/dens'.format(ineq+1)].create_group('{:05}'.format(i))
    adgafile['ineq-{:03}/magn'.format(ineq+1)].create_group('{:05}'.format(i))
    adgafile['ineq-{:03}/dens/{:05}/00001'.format(ineq+1,i)]=g4dens[...,i]
    adgafile['ineq-{:03}/magn/{:05}/00001'.format(ineq+1,i)]=g4magn[...,i]

adgafile.close()
segfile.close()

