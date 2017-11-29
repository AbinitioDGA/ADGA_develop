#!/usr/bin/env python
"""
This script computes:
  * the one-frequency susceptibility from the one-frequency
    two-particle Greens function as measured in w2dynamics,
  * the two-frequency connected susceptibility (threeleg vertex)
    from the two-frequency two-particle Greens function as measured in w2dynamics

It is possible to perform either one or both of the calculations mentioned above. 
The parameters are read from the respective file.
Additionally it is possible to use a third input file that contains the one-particle
Greens function. This feature is not yet tested, but may be useful, 
since usually we have the Greens function somewhere else in better quality.
"""


import numpy as np
import h5py
import sys
import optparse
#import matplotlib.pyplot as plt


# not very pythonic utility function to unravel the band-spin compound index
def index2component_general(Nbands, N, ind):
  b=np.zeros((N),dtype=np.int_)
  s=np.zeros((N),dtype=np.int_)
  bs=np.zeros((N),dtype=np.int_)
  # the proposed back conversion assumes the indices are
  # given form 0 to max-1
  ind_tmp = ind - 1
  tmp=np.zeros((N+1),dtype=np.int_)
  for i in xrange(0,N+1):
    tmp[i] = (2*Nbands)**(N-i)
  
  for i in xrange(N):
    bs[i] = ind_tmp/tmp[i+1]
    s[i] = bs[i]%2
    b[i] = (bs[i]-s[i])//2
    ind_tmp = ind_tmp - tmp[i+1]*bs[i]
  
  return tuple(bs),tuple(b),tuple(s)

# compute a compound index from orbital indices only. 
def component2index_band(Nbands, N, b):
  ind = 1
  for i in xrange(N):
     ind = ind + Nbands**(N-i-1)*b[i]
  return ind

def componentBS2index_general(Nbands, N, bs):

  ind = 1
  for i in xrange(N):
     ind = ind + (2*Nbands)**(N-i-1)*bs[i]
  return ind


def read_parameters(conf):
  #parameters
  if conf['do_p2']:
    h5file=conf['input_p2']
  elif conf['do_p3']:
    h5file=conf['input_p3']
  conf['beta'] = h5file[".config"].attrs["general.beta"]
  conf['nbands'] = h5file[".config"].attrs["atoms.{}.nd".format(conf['ineq'])]
  conf['niw'] = h5file[".config"].attrs["qmc.niw"]
  if conf['do_p2']:
    conf['n2iwb'] = conf['input_p2'][".config"].attrs["qmc.n2iwb"]
  if conf['do_p3']:
    conf['n3iwb'] = conf['input_p3'][".config"].attrs["qmc.n3iwb"]
    conf['n3iwf'] = conf['input_p3'][".config"].attrs["qmc.n3iwf"]
  

def read_giw(conf,options):
  nbands=conf['nbands']

  def extract_giw_worm(grp):
    giw_arr = grp['giw-worm/value'].value
    giw_arr = giw_arr.reshape((2*nbands,2*nbands,giw_arr.shape[-1]))
    niw = giw_arr.shape[-1]//2
    giw = np.zeros(shape=(2*nbands,2*niw),dtype=complex)
    #extracting diagonal from worm quantity
    for ibs in xrange(2*nbands):
      giw[ibs,:] = giw_arr[ibs,ibs,:]
    return giw

  if conf['ext_giw']:
    ineq_grp=conf['input_giw'][options.giw_iter+('/ineq-{:03}'.format(conf['ineq']))]
    keys=ineq_grp.keys()
    if 'giw-worm' in keys:
      giw=extract_giw_worm(ineq_grp)
    elif 'giw' in keys:
      giw=ineq_grp['giw/value'].value.reshape((2*nbands,-1))
    else:
      sys.exit('Error: no giw found in giw file!')

  elif conf['do_p2']:
    ineq_grp=conf['input_p2']['stat-001/ineq-{:03}'.format(conf['ineq'])]
    keys=ineq_grp.keys()
    if 'giw-worm' in keys:
      giw=extract_giw_worm(grp)
    elif 'giw' in keys:
      giw=ineq_grp['giw/value'].value.reshape((2*nbands,-1))

  elif conf['do_p3']:
    ineq_grp=conf['input_p3']['stat-001/ineq-{:03}'.format(conf['ineq'])]
    keys=ineq_grp.keys()
    if 'giw-worm' in keys:
      giw=extract_giw_worm(grp)
    elif 'giw' in keys:
      giw=ineq_grp['giw/value'].value.reshape((2*nbands,-1))

  return giw 


def get_gg_ph(giw,niw=None,n3iwf=None,n3iwb=None,nbands=None,**kwargs):
  gg=np.zeros((2*nbands,2*nbands,2*n3iwf,2*n3iwb+1),dtype=np.complex)
  
  for iwb in xrange(-n3iwb,n3iwb+1):
    gg[...,iwb+n3iwb]=giw[:,None,niw-n3iwf:niw+n3iwf]*giw[None,:,niw-n3iwf-iwb:niw+n3iwf-iwb]

  return gg


def slice_center(subsize, size):
    """Returns a slice for the center elements of a bigger array"""
    if subsize < 0 or size < 0:
        raise ValueError("sizes must be bigger than zero")
    if subsize > size: 
        raise ValueError("subsize > size")
    if subsize % 2 != size % 2:
        raise ValueError("size and subsize are different MOD 2")
    if subsize == size:
        return slice(None)
    else:
        startidx = (size - subsize)//2
        cslice = slice(startidx, startidx + subsize)
        return cslice


def read_p3(ineq=None,input_p3=None,nbands=None,**kwargs):
  if 'stat-001' in input_p3.keys():
    p3 = input_p3["stat-001/ineq-{:03}/p3iw-worm/value".format(ineq)].value
  elif 'worm-001' in input_p3.keys():
    base_gr=input_p3["worm-001/ineq-{:03}/p3iw-worm".format(ineq)]
    groups=base_gr.keys()
    nf,nb=base_gr[groups[0]+'/value'].shape
    nf_small,nb_small=nf//2,nb//2+1
    #p3=np.zeros((2*nbands,2*nbands,2*nbands,2*nbands,nf_small,nb_small),dtype=np.complex) # nf,nb
    p3=np.zeros((2*nbands,2*nbands,2*nbands,2*nbands,nf,nb),dtype=np.complex) # nf,nb
    for gr in groups:
      bs,_,_=index2component_general(nbands,4,int(gr))
      #p3[bs]=base_gr[gr+'/value'].value[nf//2-nf_small//2:nf//2+nf_small//2,nb//2-nb_small//2:nb//2+nb_small//2+1]
      p3[bs]=base_gr[gr+'/value'].value
    return p3



parser = optparse.OptionParser(usage = "%prog [OPTIONS] <hdf-file>")
parser.add_option("--p2-file",dest="p2_file",help="w2dynamics hdf5 file containing one-frequency two-particle greens function")
parser.add_option("--p3-file",dest="p3_file",help="w2dynamics hdf5 file containing two-frequency two-particle greens function")
parser.add_option("--giw-file",dest="giw_file",help="w2dynamics hdf5 file containing one-particle greens function")
parser.add_option("--giw-iter",dest="giw_iter",help="iteration name, where giw can be found in external file, e.g. stat-001")
parser.add_option("--threeleg-file",dest="threeleg_file",help="output file name for threeleg vertex")
parser.add_option("--susc-file",dest="susc_file",help="output file name for one-frequency susceptibility")
parser.add_option("--ineq",dest="ineq",default=1,type=int,help="number of inequivalent atom, larger or equal to 1.")

options, args = parser.parse_args()



conf={}
conf['do_p2']=False
conf['do_p3']=False
conf['ineq'],ineq=options.ineq,options.ineq

if options.p2_file is not None:
  conf['input_p2']=h5py.File(options.p2_file,'r')
  conf['do_p2']=True
if options.p3_file is not None:
  conf['input_p3']=h5py.File(options.p3_file,'r')
  conf['do_p3']=True
if options.giw_file is not None:
  conf['ext_giw']=True
  conf['input_giw']=h5py.File(options.giw_file,'r')
else:
  conf['ext_giw']=False


read_parameters(conf)
nbands=conf['nbands']

giw=read_giw(conf,options)
conf['niw']=giw.shape[-1]//2
niw=conf['niw']

if conf['do_p2']:# subtract density terms to get susceptibility
  print 'processing p2...'
  p2 = conf['input_p2']["stat-001/ineq-{:03}/p2iw-worm/value".format(ineq)].value
  n2iwb=p2.shape[0]//2
  conf['n2iwb']=n2iwb
  # density 
  occ = conf['input_p2']["stat-001/ineq-{:03}/occ/value".format(ineq)].value
  dens = np.diagonal(np.diagonal(occ,0,-4,-2),0,-3,-2)
  dens = dens.reshape(2*nbands)
  k1 = p2
  # kernel 1 functions - subtract density terms to obtain susceptibility
  for i in xrange(2*nbands):
     for k in xrange(2*nbands):
        k1[i,i,k,k,n2iwb] -= conf['beta']*(1-dens[i])*(1-dens[k])



  # we save the one-frequency susceptibilities for convenience
  f_susc=h5py.File(options.susc_file,'w')

  for bs1 in xrange(2*nbands):
    for bs2 in xrange(2*nbands):
      for bs3 in xrange(2*nbands):
        for bs4 in xrange(2*nbands):
          if np.any(k1[bs1,bs2,bs3,bs4]!=0.):
             gr_name='ineq-{:03}/{:05}'.format(ineq,componentBS2index_general(nbands,4,[bs1,bs2,bs3,bs4]))
             f_susc[gr_name]=k1[bs1,bs2,bs3,bs4]

  f_susc.close()

    

if conf['do_p3']:# subtract all disconnected terms to get threeleg vertex
  print 'processing p3...'
  #p3 = conf['input_p3']["worm-001/ineq-001/p3iw-worm/value"].value
  p3 = read_p3(**conf)
  n3iwf,n3iwb=p3.shape[-2:]
  n3iwf=n3iwf//2
  n3iwb=n3iwb//2
  conf['n3iwf']=n3iwf
  conf['n3iwb']=n3iwb
  gg_ph=get_gg_ph(giw,**conf)
  # density 
  try:
    occ = conf['input_p3']["stat-001/ineq-{:03}/occ/value".format(ineq)].value
  except KeyError:
    occ = conf['input_giw']["dmft-last/ineq-{:03}/occ/value".format(ineq)].value
  dens = np.diagonal(np.diagonal(occ,0,-4,-2),0,-3,-2)
  dens = dens.reshape(2*nbands)
  print dens

  k2 = np.zeros_like(p3,dtype=complex)
  k2 = -conf['beta']*p3

  # subtract disconnected terms
  for i in xrange(2*nbands):
    for j in xrange(2*nbands):
      for k in xrange(2*nbands):
        for l in xrange(2*nbands):
              
          # straight terms
          if i==j and k==l:
            k2[i,j,k,l,:,n3iwb] -= giw[i,slice_center(2*n3iwf,2*niw)]*(1-dens[k])*conf['beta']
  
          # cross terms
          if i==l and j==k:
            k2[i,j,k,l,:,:] -= gg_ph[i,k]


  # amputate 2 legs
  for i in xrange(2*nbands):
    for j in xrange(2*nbands):
      for k in xrange(2*nbands):
        for l in xrange(2*nbands):
          k2[i,j,k,l,...]=np.divide(k2[i,j,k,l,...],gg_ph[i,j])
           
  
  fout=h5py.File(options.threeleg_file,'a')
  for bs1 in xrange(2*nbands):
    for bs2 in xrange(2*nbands):
      for bs3 in xrange(2*nbands):
        for bs4 in xrange(2*nbands):
          if np.any(k2[bs1,bs2,bs3,bs4]!=0.):
             gr_name='ineq-{:03}/{:05}'.format(ineq,componentBS2index_general(nbands,4,[bs1,bs2,bs3,bs4]))
             fout[gr_name]=k2[bs1,bs2,bs3,bs4]
  fout.close()
