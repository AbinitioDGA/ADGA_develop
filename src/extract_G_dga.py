#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import h5py
import sys

f = h5py.File('Output-insulating-u50-transposed/adga-20171113-184415.568-output.hdf5','r')

ndim = f['input/hk'].shape[0]
beta = f['input/beta'][()]
nqx, nqy, nqz = f['input/nqpxyz'][0], f['input/nqpxyz'][1], f['input/nqpxyz'][2]
iwf = f['input/iwfmax_small'][()]
iwfdmft = f['input/siw'].shape[-1]

fmats = np.linspace(-(iwf*2-1)*np.pi/beta,(iwf*2-1)*np.pi/beta,2*iwf) # fermionic matsubara axis

# load h5py datasets into predefined numpy arrays
siwdga = np.zeros((ndim,ndim,nqx,nqy,nqz,2*iwf),dtype=np.complex128)
f['selfenergy/nonloc/dga'].read_direct(siwdga)
siwdgaloc = np.zeros((ndim,ndim,2*iwf),dtype=np.complex128)
f['selfenergy/loc/dga_ksum'].read_direct(siwdgaloc)
siw = np.zeros((ndim,iwfdmft),dtype=np.complex128)
f['input/siw'].read_direct(siw)
giw = np.zeros((ndim,iwfdmft),dtype=np.complex128)
f['input/giw'].read_direct(giw)

# reading done
f.close()


# local DGA selfenergy
g0iwinv = np.zeros_like(giw, dtype=np.complex128)
g0iwinv = giw[...]**(-1) + siw[...]

# G0^(-1) = G_dmft^(-1) + Sigma_dmft^(-1) -> orbital diagonal
# G_dga^(-1)(k) = G0^(-1) - Sigma_dga(k)
# G_dga(k) = (G_dga^(-1)(k))^(-1)
# sum over k-points
# plot diagonal elements of that
if True:
    gdgainv = np.zeros_like(siwdga,dtype=np.complex128)
    gdga = np.zeros_like(siwdga,dtype=np.complex128)

    gdgainv = -siwdga

    for iband in xrange(ndim):
        gdgainv[iband,iband,...] += g0iwinv[iband,None,None,None,iwfdmft/2-iwf:iwfdmft/2+iwf]

    for kx in xrange(nqx):
        for ky in xrange(nqy):
            for kz in xrange(nqz):
                print(kx,ky,kz)
                for iw in xrange(2*iwf):
                    gdga[:,:,kx,ky,kz,iw] = scipy.linalg.inv(gdgainv[:,:,kx,ky,kz,iw])

    gdgaksum = np.zeros((ndim,ndim,2*iwf),dtype=np.complex128)
    gdgaksum = np.sum(gdga,axis=(2,3,4))/(nqx*nqy*nqz) # divide by number of k-points

    for iband in xrange(ndim):
        plt.plot(fmats,giw[iband,iwfdmft/2-iwf:iwfdmft/2+iwf].imag,'-x', label='dmftband: {}'.format(iband+1))
        plt.plot(fmats,gdgaksum[iband,iband,:].imag,'--o',label='dgaband: {}'.format(iband+1))

    plt.legend()
    plt.show()
    raw_input()

# G0^(-1) = G_dmft^(-1) + Sigma_dmft^(-1) -> orbital diagonal
# G_dga^(-1) = G0^(-1) - Sigma_dga_loc
# G_dga = (G_dga^(-1))^(-1)
# plot diagonal elements of that
if False:
    gdgainv = np.zeros_like(siwdgaloc)
    gdgainv = -siwdgaloc
    for iband in xrange(ndim):
        gdgainv[iband,iband,:] += g0iwinv[iband,iwfdmft/2-iwf:iwfdmft/2+iwf]

    gdga = np.zeros((ndim,ndim,2*iwf), dtype=np.complex128)

    for iw in xrange(2*iwf):
        gdga[...,iw] = scipy.linalg.inv(gdgainv[...,iw])

    for iband in xrange(ndim):
        plt.plot(fmats,giw[iband,iwfdmft/2-iwf:iwfdmft/2+iwf].imag, '-x', label='dmftband: {}'.format(iband+1))
        plt.plot(fmats,gdga[iband,iband,:].imag,'--o',label='dgaband: {}'.format(iband+1))

    plt.legend()
    plt.show()
    raw_input()
