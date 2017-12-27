#!/usr/bin/env python

import h5py
import numpy as np
import sys

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

  return bs,b,s

# compute a compound index from orbital indices only.
def component2index_band(Nbands, N, b):
  ind = 1
  for i in xrange(N):
     ind = ind + Nbands**(N-i-1)*b[i]
  return ind

def ask_for_input():
  conf={}
  print "Which quantity do you want to symmetrize?"
  targ=int(raw_input("1 -> 1 frequency (b), 2 -> 2 frequencies (fb), 3 -> 3 frequencies (ffb): "))
  if targ==1:
    conf['target']='1freq'
  elif targ==2:
    conf['target']='2freq'
  elif targ==3:
    conf['target']='3freq'

  nineq=int(raw_input('Number of inequivalent atoms: '))
  conf['nineq']=nineq
  filename_nosym=raw_input('Filename of the not symmetrized data: ')
  conf['infile']=filename_nosym
  Nbands=map(int,raw_input('Number of correlated orbitals in each inequivalent atom (seperated by spaces): ').split())
  conf['Nbands']=Nbands
  if len(Nbands) != nineq:
    print 'Provide number of orbitals for each inequivalent atom!'
    sys.exit()
  filename_sym=raw_input('Outputfile for symmetrized Vertex: ')
  conf['outfile']=filename_sym
  conf['sym_type']=raw_input('SU2 symmetry only (s) or SU2 AND orbital symmetry (o)?: ')
  return conf

def get_groups(infile='infile.hdf5',Nbands=[1],target=1,nineq=1,**kwargs):
  groups=[]
  bgroups=[]
  for ineq in xrange(nineq):
    f=h5py.File(infile,'r')
    if target=='1freq':
      gr_str=f['ineq-{:03}'.format(ineq+1)].keys()
    elif target=='2freq':
      gr_str=f['ineq-{:03}'.format(ineq+1)].keys()
    elif target=='3freq':
      gr_str=f['worm-last/ineq-{:03}/g4iw-worm'.format(ineq+1)].keys()
    f.close()
    groups.append([])
    bgroups.append([])
    groups[ineq]=[]
    bgroups[ineq]=[]
    bind=[]
    for gr in gr_str:
      bs,b,s=index2component_general(Nbands[ineq],4,int(gr))
      groups[ineq].append({'group':int(gr),'spins':tuple(s),'bands':tuple(b),'band-spin':tuple(bs)})
      bgr=component2index_band(Nbands[ineq],4,b)
      if not bgr in bind:
        bgroups[ineq].append({'bgroup':bgr,'bands':tuple(b)})
      bind.append(bgr)
  conf['groups']=groups
  conf['bgroups']=bgroups

def get_fbox(infile=None,target=None,**kwargs):
  f=h5py.File(infile,'r')
  if target=='1freq':
    n2iwb=f['ineq-001/00001'].shape[0]
    print n2iwb
    conf['n2iwb']=n2iwb//2
  elif target=='2freq':
    n3iwf,n3iwb = f['ineq-001/00001'].shape
    conf['n3iwf'],conf['n3iwb'] = n3iwf//2,n3iwb//2
    n3iwf,n3iwb=conf['n3iwf'],conf['n3iwb']
  elif target=='3freq':
    # new worm format -- f f b
    _,n4iwf,n4iwb = f['worm-last/ineq-001/g4iw-worm/00001/value'].shape # always here
    conf['n4iwf']=n4iwf//2
    conf['n4iwb']=n4iwb//2
    n4iwf,n4iwb=conf['n4iwf'],conf['n4iwb']

  f.close()


def initialize_output(f1,h5f,bgroups=None,nineq=None,n2iwb=None,n3iwf=None,n3iwb=None,n4iwf=None,n4iwb=None,target=None,**kwargs):
  for ineq in xrange(nineq):
    dset_ineq=h5f.create_group('ineq-{:03}'.format(ineq+1))
    dset_dens=dset_ineq.create_group('dens')
    dset_magn=dset_ineq.create_group('magn')

    if target=='1freq':
      for bgr in [d['bgroup'] for d in bgroups[ineq]]:
        dset_dens['{:05}'.format(bgr)]=np.zeros((2*n2iwb+1),dtype=np.complex128)
        dset_magn['{:05}'.format(bgr)]=np.zeros((2*n2iwb+1),dtype=np.complex128)
    elif target=='2freq':
      for iwb in xrange(2*n3iwb+1):
        dset_d1=dset_dens.create_group('{:05}'.format(iwb))
        dset_m1=dset_magn.create_group('{:05}'.format(iwb))
        for bgr in [d['bgroup'] for d in bgroups[ineq]]:
          dset_d1['{:05}'.format(bgr)]=np.zeros((2*n3iwf),dtype=np.complex128)
          dset_m1['{:05}'.format(bgr)]=np.zeros((2*n3iwf),dtype=np.complex128)
    elif target=='3freq':
      if ineq==0:
        f1.copy('.axes',h5f)
      for iwb in xrange(2*n4iwb+1):
        dset_d1=dset_dens.create_group('{:05}'.format(iwb))
        dset_m1=dset_magn.create_group('{:05}'.format(iwb))
        for bgr in [d['bgroup'] for d in bgroups[ineq]]:
          dset_d1['{:05}/value'.format(bgr)]=np.zeros((2*n4iwf,2*n4iwf),dtype=np.complex128)
          dset_m1['{:05}/value'.format(bgr)]=np.zeros((2*n4iwf,2*n4iwf),dtype=np.complex128)


def get_symgroups(ch,gr,nd,sym_type=None,**kwargs):
  if ch == 'dens':
    if gr['spins'] in [(0,0,0,0),(0,0,1,1),(1,1,0,0),(1,1,1,1)]:
      action = '+'
    elif gr['spins'] in [(0,1,1,0),(1,0,0,1)]:
      action = '0'
    else:
      print 'unknown spin combination'
      sys.exit()
  elif ch == 'magn':
    if gr['spins'] in [(0,0,0,0),(0,1,1,0),(1,0,0,1),(1,1,1,1)]:
      action = '+'
    elif gr['spins'] in [(0,0,1,1),(1,1,0,0)]:
      action = '-'
    else:
      print 'unknown spin combination'
      sys.exit()

  if sym_type=='o':
    b1,b2,b3,b4=gr['bands']
    symgroups=[]
    if b1==b2 and b1==b3 and b1==b4:
      for i in xrange(nd):
        symgroups.append(component2index_band(nd,4,[i,i,i,i]))
    elif b1==b2 and b3==b4 and b2!=b3:
      for i in xrange(nd):
        for j in xrange(nd):
          if i!=j:
            symgroups.append(component2index_band(nd,4,[i,i,j,j]))
    elif b1==b3 and b2==b4 and b1!=b2:
      for i in xrange(nd):
        for j in xrange(nd):
          if i!=j:
            symgroups.append(component2index_band(nd,4,[i,j,i,j]))
    elif b1==b4 and b2==b3 and b1!=b2:
      for i in xrange(nd):
        for j in xrange(nd):
          if i!=j:
            symgroups.append(component2index_band(nd,4,[i,j,j,i]))

  else:
    symgroups=[component2index_band(nd,4,gr['bands'])]

  return action,symgroups

def read_and_add(h5in,h5out,ineq,igr,channel,action,symgroups,target=None,n3iwb=None,n4iwb=None,**kwargs):
  if channel=='dens':
    prefactor = 2.0
  elif channel=='magn':
    prefactor = 4.0
  else:
    print 'unkown channel'
    sys.exit()

  if target=='1freq':
    x = h5in['ineq-{:03}/{:05}'.format(ineq+1,igr)].value/float(prefactor*len(symgroups))
    for gr in symgroups:
      if action=='+':
        h5out['ineq-{:03}/{}/{:05}'.format(ineq+1,channel,gr)][...]+=x
      elif action=='-':
        h5out['ineq-{:03}/{}/{:05}'.format(ineq+1,channel,gr)][...]-=x
      elif action=='0':
        pass
  elif target=='2freq':
    x = h5in['ineq-{:03}/{:05}'.format(ineq+1,igr)].value/float(prefactor*len(symgroups))
    for iwb in xrange(2*n3iwb+1):
      for gr in symgroups:
        if action=='+':
          h5out['ineq-{:03}/{}/{:05}/{:05}'.format(ineq+1,channel,iwb,gr)][...]+=x[:,iwb]
        elif action=='-':
          h5out['ineq-{:03}/{}/{:05}/{:05}'.format(ineq+1,channel,iwb,gr)][...]-=x[:,iwb]
        elif action=='0':
          pass
  elif target=='3freq':
    x = h5in['worm-last/ineq-{:03}/g4iw-worm/{:05}/value'.format(ineq+1,igr)].value/float(prefactor*len(symgroups))
    for iwb in xrange(2*n4iwb+1):
      for gr in symgroups:
        if action=='+':
          h5out['ineq-{:03}/{}/{:05}/{:05}/value'.format(ineq+1,channel,iwb,gr)][...]+=x[...,iwb].transpose()
        elif action=='-':
          h5out['ineq-{:03}/{}/{:05}/{:05}/value'.format(ineq+1,channel,iwb,gr)][...]-=x[...,iwb].transpose()
        elif action=='0':
          pass



#===============================================================================
#================================ Script start =================================
#===============================================================================


conf=ask_for_input()
#conf={'nineq': 1, 'target': '3freq', 'sym_type': 'o', 'outfile': 'out.hdf5', 'Nbands': [3,3], 'infile': 'vertex_full_newformat.hdf5'}
print conf


get_groups(**conf)
get_fbox(**conf)

f1=h5py.File(conf['infile'],'r')
f2=h5py.File(conf['outfile'],'w')

initialize_output(f1,f2,**conf)

for ineq in xrange(conf['nineq']):
  for ch in ['dens','magn']:
      for gr in conf['groups'][ineq]:
        action,symgroups=get_symgroups(ch,gr,conf['Nbands'][ineq],**conf)
        print 'group {},'.format(gr['group']),'channel: {},'.format(ch),'action: {},'.format(action),'{} equivaluent band groups:'.format(len(symgroups)),symgroups
        read_and_add(f1,f2,ineq,gr['group'],ch,action,symgroups,**conf)

f1.close()
f2.close()
