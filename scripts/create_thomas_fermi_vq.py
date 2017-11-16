#! /usr/bin/env python

from __future__ import print_function, division
import h5py
import numpy as np
import matplotlib.pyplot as plt

k_tf = 0.3 # 0.3 #2.5 # k thomas fermi
nqx, nqy, nqz = 20, 20, 20 # q points
scx, scy, scz =  2*np.pi, 2*np.pi, 2*np.pi # actual lengths of BZ

def q_value(iq):
    ''' return q value for q-point
    iq = iqz + iqy*nqz + iqx*nqy*nqz '''
    iqx=iq//(nqy*nqz)
    iqy=(iq-iqx*nqy*nqz)//nqz
    iqz=(iq-iqy*nqz-iqx*nqy*nqz)
    # scale up so BZ is defined from [0, 2pi)
    return scx*iqx/nqx, scy*iqy/nqy, scz*iqz/nqz

def q_point(iq):
    ''' return q-point
    iq = iqz + iqy*nqz + iqx*nqy*nqz '''
    iqx=iq//(nqy*nqz)
    iqy=(iq-iqx*nqy*nqz)//nqz
    iqz=(iq-iqy*nqz-iqx*nqy*nqz)
    return iqx,iqy,iqz

def tf_value(iqx,iqy,iqz):
    return 4*np.pi / (iqx**2+iqy**2+iqz**2+k_tf**2)

dset_axes = np.zeros((nqx*nqy*nqz, 3), dtype=np.float64)
dset_vq = np.zeros((nqx*nqy*nqz, ), dtype=np.complex128)

# this should be periodical in the BZ ...
# should mean that
# charge at q = 0 0 0, q = 2pi 0 0 , q = 0 2pi 0, 2pi 2pi 0
# charge at q = 0 0 2pi, q = 2pi 0 2pi, q = 0 2pi 2pi, 2pi 2pi 2pi

# if k_tf is small -> we need more neighbours
# 8 corners 3 directions per corner to maintain symmetry
# -2pi 0 0, 0 -2pi 0, 0 0 -2pi
# 4pi 0 0, 2pi -2pi 0, 2pi 0 -2pi
# 0 4pi 0, -2pi 2pi 0, 0 2pi -2pi
# 4pi 2pi 0, 2pi 4pi 0, 2pi 2pi -2pi
# -2pi 0 2pi, 0 -2pi 2pi, 0 0 4pi
# 4pi 0 2pi, 2pi -2pi 2pi, 2pi 0 4pi
# 0 4pi 2pi, -2pi 2pi 2pi, 0 2pi 4pi
# 4pi 2pi 2pi, 2pi 4pi 2pi, 2pi 2pi 4pi

for iq in xrange(nqx*nqy*nqz):
    iqx, iqy, iqz = q_value(iq)
    dset_axes[iq,:] = np.array([iqx, iqy, iqz], dtype=np.float64)
    # first 8 : 8 corners of BZ
    # next 24 : 3 * 8 neighbours of those
    dset_vq[iq] = tf_value(iqx,iqy,iqz) \
                + tf_value(scx-iqx,iqy,iqz) \
                + tf_value(iqx,scy-iqy,iqz) \
                + tf_value(scx-iqx,scy-iqy,iqz) \
                + tf_value(iqx,iqy,scz-iqz) \
                + tf_value(scx-iqx,iqy,scz-iqz) \
                + tf_value(iqx,scy-iqy,scz-iqz) \
                + tf_value(scx-iqx,scy-iqy,scz-iqz) \
                + tf_value(scx+iqx,iqy,iqz) \
                + tf_value(iqx,scy+iqy,iqz) \
                + tf_value(iqx,iqy,scz+iqz) \
                + tf_value(2*scx-iqx,iqy,iqz) \
                + tf_value(scx-iqx,scy+iqy,iqz) \
                + tf_value(scx-iqx,iqy,scz+iqz) \
                + tf_value(iqx,2*scy-iqy,iqz) \
                + tf_value(scx+iqx,scy-iqy,iqz) \
                + tf_value(iqx,scy-iqy,scz+iqz) \
                + tf_value(2*scx-iqx,scy-iqy,iqz) \
                + tf_value(scx-iqx,2*scy-iqy,iqz) \
                + tf_value(scx-iqx,scy-iqy,scz+iqz) \
                + tf_value(scx+iqx,iqy,scz-iqz) \
                + tf_value(iqx,scy+iqy,scz-iqz) \
                + tf_value(iqx,iqy,2*scz-iqz) \
                + tf_value(2*scx-iqx,iqy,scz-iqz) \
                + tf_value(scx-iqx,scy+iqy,scz-iqz) \
                + tf_value(scx-iqx,iqy,2*scz-iqz) \
                + tf_value(iqx,2*scy-iqy,scz-iqz) \
                + tf_value(scx+iqx,scy-iqy,scz-iqz) \
                + tf_value(iqx,scy-iqy,2*scz-iqz) \
                + tf_value(2*scx-iqx,scy-iqy,scz-iqz) \
                + tf_value(scx-iqx,2*scy-iqy,scz-iqz) \
                + tf_value(scx-iqx,scy-iqy,2*scz-iqz)

# subtract mean so vq is purely nonlocal
dset_vq = dset_vq - np.mean(dset_vq)

# h5py stuff
f=h5py.File('Vq_hubbard_tf.hdf5','w')
axes=f.create_group('.axes')

# for the mean time only with 1 band -> 1 1 1 1 -> 00001
axes.create_dataset('Q-points',data=dset_axes)
f.create_dataset('00001',data=dset_vq)

print('maximum with offset of U=2: ', np.max(dset_vq.real+2))
print('minimum with offset of U=2: ', np.min(dset_vq.real+2))
plt.plot(dset_vq.real+2)
plt.plot(np.full_like(dset_vq,2,dtype=np.float64))
plt.plot(np.full_like(dset_vq,0,dtype=np.float64))
plt.show()


# plotting extended stuff

dset_cut = np.zeros((nqx,nqy,nqz), dtype=np.float64)
for iq in xrange(nqx*nqy*nqz):
    iqx, iqy, iqz = q_point(iq)
    dset_cut[iqx,iqy,iqz] = dset_vq[iq].real

# test for plotting (extend to full BZ)
full_bz=np.zeros((nqx+1,nqy+1,nqz+1),dtype=np.float64)
full_bz[:-1,:-1,:-1]=dset_cut[...]#-dset2_cut[...]
full_bz[:,:,nqz]=full_bz[:,:,0]
full_bz[:,nqy,:]=full_bz[:,0,:]
full_bz[nqx,:,:]=full_bz[0,:,:]
full_bz[nqx,nqy,:]=full_bz[0,0,:]
full_bz[nqx,:,nqz]=full_bz[0,:,0]
full_bz[:,nqy,nqz]=full_bz[:,0,0]
full_bz[nqx,nqy,nqz]=full_bz[0,0,0]

plt.pcolormesh(full_bz[:,:,0],cmap='seismic')
plt.tight_layout()
plt.colorbar()
plt.show()
