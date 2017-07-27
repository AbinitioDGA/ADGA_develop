#!/usr/bin/env python
'''
This script reads one HDF5 file that contains nonlocal interaction U(r) and a list of r vectors.

It writes one file with extrapolated U(r) and one file with Fouriertransformed V(q).
The constant U(0) is currently not subtracted. 

It is necessary to provide one input and two output filenames,
and the number of r points for the extrapolated U(r) in x-, y- and z-direction,
and the number of q points for the Fourier transform in qx-, qy- and qz-direction.

The type of the fit function cannot be specified interactively, but only in the script,
approx. line 100:
                 fitfun=[yukawa | simple_exp]


sample execution line:
./V_r_extrapol.py  ~/Data/SVO/U-00RR-vec U_r_extrapol.dat V_r.hdf5 10 10 10 10 10 10



Author: Josef Kaufmann,
last change: May 2017
'''

import numpy as np
from scipy.optimize import curve_fit



import sys
import optparse



__version__ = "1.0"


parser = optparse.OptionParser(usage = "%prog [OPTIONS] input output nx ny nz",
                               version = ".".join(map(str,__version__)))

#parser.add_option("--ncomp", dest="ncomp",default=0,
#                 help="number of U-matrix components contained in the input file. (i.e. number of columns minus 3)")

parser.add_option("--nd", dest="nd",default=0, help="number of orbitals")
parser.add_option("--interactive",dest="show_plots",default=False,action="store_true", help="switch on interactive plots")


options, args = parser.parse_args()
if len(args) != 9:
  parser.error("expected input and output filename and number of r and q points in each direction ")


if options.show_plots:
  from matplotlib import pyplot as plt


filename_in=args[0]
filename_out=args[1]
outfilename_vq=args[2]
nx=int(args[3])
ny=int(args[4])
nz=int(args[5])
nqx=int(args[6])
nqy=int(args[7])
nqz=int(args[8])
nr=(2*nx+1)*(2*ny+1)*(2*nz+1)
#ncomp=int(options.ncomp) # number of components of U(r)
nd=int(options.nd) # number of orbitals just needed to write output file.

#print nd,'Bands'
#print 'U matrix has',ncomp,'components'


# read the input data from text file
#infile=open(filename_in,'r')
#data=infile.readlines()
#n = len(data)
#r_vec=np.zeros((n,3),dtype=np.float_)
#r_abs=np.zeros((n),dtype=np.float_)
#U_r_in=np.zeros((n,ncomp),dtype=np.float_)
#for i in xrange(n):
#  r_vec[i]=data[i].split()[:3]
#  r_abs[i]=np.linalg.norm(r_vec[i])
#  U_r_in[i]=np.float_(data[i].split()[-ncomp:])
#infile.close()

# read the input data from HDF5 file
import h5py
infile=h5py.File(filename_in,'r')
groups=[key for key in infile.keys() if key.isdigit()]
ncomp=len(groups)
print 'Nonlocal interaction has',ncomp,'components'
r_vec=infile['.axes/R-points'].value
r_abs=np.array(map(np.linalg.norm,r_vec))
print r_abs.shape
n=r_vec.shape[0]
U_r_in=np.zeros((n,ncomp),dtype=np.float)
for igr,gr in enumerate(groups):
  U_r_in[:,igr]=infile[gr+'/value'].value


# plot the input data vs abs(r)
plotcomp=0 # component to plot
if options.show_plots:
  f1=plt.figure(1)
  plt.plot(r_abs,U_r_in[:,plotcomp].real,marker='.',linestyle='None')
  f1.show()


# do the fit
# since the function we use cannot do a constrained optimization,
# I add a "penalization" if b<=0.
# I very much hope that this is correct.
def yukawa_2(r,a,b):
  penalization = 0.
  if b<=0.:
    penalization=1.e20
  return a*np.exp(-b*np.abs(r))/np.abs(r) + penalization

def yukawa_3(r,a,b,c):
  penalization = 0.
  if b<=0.:
    penalization=1.e20
  return a*np.exp(-b*np.abs(r))/np.abs(r) + penalization + c


def simple_exp(r,a,b):
  penalization = 0.
  if b<=0.:
    penalization=1.e20
  return a*np.exp(-b*np.abs(r)) + penalization
#  return 1.95*np.exp(-b*np.abs(r))


fitfun=yukawa_2

if fitfun==yukawa_2 or fitfun==yukawa_3:
# if we fit to yukawa, take out the value at r=0, because we cannot fit it otherwise.
  ydata=np.zeros((n-1,ncomp),dtype=np.float_)
  rdata=np.zeros((n-1),dtype=np.float_)
  j=0
  for i in xrange(n):
    if r_abs[i] >= 0.001:
    #if r_abs[i] >=3.:
      ydata[j]=U_r_in[i]
      rdata[j]=r_abs[i]
      j=j+1
  ydata=ydata[:j-1]
  rdata=rdata[:j-1]
else:
  ydata=np.zeros((n,ncomp),dtype=np.float_)
  rdata=np.zeros((n),dtype=np.float_)
  j=0
  for i in xrange(n):
    ydata[j]=U_r_in[i]
    rdata[j]=r_abs[i]
    j=j+1

n_fitparam=2
if fitfun==yukawa_3:
  n_fitparam=3
fitparam=np.zeros((ncomp,n_fitparam),dtype=float)
fitcov = np.zeros((ncomp,n_fitparam,n_fitparam),dtype=float)
for iu in xrange(ncomp):
  fitparam[iu], fitcov[iu] = curve_fit(fitfun,rdata,ydata[:,iu])
  print 'component',iu
  print 'opt. fit params.',fitparam[iu]
  print 'covariance',fitcov[iu]
  if n_fitparam==3:
    fitparam[iu,2]=0.

if options.show_plots:
  f2=plt.figure(2)
  plt.plot(np.linspace(min(rdata),max(rdata)),fitfun(np.linspace(min(rdata),max(rdata)),*fitparam[plotcomp]))
  plt.plot(r_abs,U_r_in[:,plotcomp].real,marker='.',linestyle='None')
  plt.xlabel('r')
  plt.ylabel('V(r)')
  plt.savefig('V_fit.png',dpi=300)
  f2.show()
print fitparam[plotcomp],min(rdata),max(rdata),fitfun(max(rdata),*fitparam[plotcomp])
#for iplt in xrange(ncomp):
#  plt.figure(iplt+10)
#  plt.plot(np.linspace(min(rdata),max(rdata)),fitfun(np.linspace(min(rdata),max(rdata)),fitparam[iplt,0],fitparam[iplt,1]))
#  plt.plot(r_abs,U_r_in[:,iplt].real,marker='.',linestyle='None')
#  plt.show()
#raw_input()
#sys.exit()
# Now we evaluate the fit function on a grid specified by the command line arguments

# First we need the lattice spacing
xvals=np.unique(r_vec[:,0])
yvals=np.unique(r_vec[:,1])
zvals=np.unique(r_vec[:,2])
a=xvals[1]-xvals[0]
b=yvals[1]-yvals[0]
c=zvals[1]-zvals[0]

print 'lattice spacing ',a,b,c

r_ext=np.zeros((2*nx+1,2*ny+1,2*nz+1,3),dtype=np.float_)
U_ext=np.zeros((2*nx+1,2*ny+1,2*nz+1,ncomp),dtype=np.float_)


outfile=open(filename_out,'w')
outfile.write('{} {} {} {} {} \n'.format(nr,nd,a,b,c))
for i in xrange(2*nx+1):
  print i
  for j in xrange(2*ny+1):
    for k in xrange(2*nz+1):
      r_ext[i,j,k,:]=[(i-nx)*a,(j-ny)*b,(k-nz)*c]
      for iu in xrange(ncomp):
        U_ext[i,j,k,iu]=fitfun(np.linalg.norm(r_ext[i,j,k]),*fitparam[iu])
      if np.linalg.norm(r_ext[i,j,k,:]) < np.max(r_abs):
        for l in xrange(n):
          if np.allclose(r_vec[l],r_ext[i,j,k],atol=np.min([a,b,c])/10) and np.linalg.norm(r_vec[l])<0.001:
            U_ext[i,j,k]=U_r_in[l]
            break
# 
      outfile.write(('{:08f}\t{:08f}\t{:08f}\n'+ncomp*'{:14f}\t'+'\n').format(r_ext[i,j,k,0],r_ext[i,j,k,1],r_ext[i,j,k,2],*U_ext[i,j,k]))
outfile.close()

if options.show_plots:
  f3=plt.figure(3)
  plt.plot(r_ext[:,ny,nz,0],U_ext[:,ny,nz,plotcomp],marker='.',linestyle='None',label='U(x)')
  # The plot in z dirction is done separately, because in HgBa2CuO4 we have different lattice spacing
  # and I want to show that the fit works also in z direction. (So it is actually unnecessary for SVO)
  #plt.plot(r_ext[nx,ny,:,2],U_ext[nx,ny,:,plotcomp],marker='*',linestyle='None',label='U(z)')
  plt.legend()
  plt.xlim(-10,10)
  plt.xlabel('x,z')
  f3.show()

U_ext=U_ext.reshape((-1,ncomp))
r_ext=r_ext.reshape((-1,3))
r_ext_h5=h5py.File('r_ext.hdf5','w')
for igr,gr in enumerate(groups):
  r_ext_h5[gr]=U_ext[...,igr]
r_ext_h5.create_group('.axes')
r_ext_h5['.axes/R-points']=r_ext



print 'Extrapolation done'
print 'Now Fourier transform...'
# TODO use a FFT routine here (this is currently not done in order to have more flexibility with the grids.)
Vq=np.zeros((nqx,nqy,nqz,ncomp),dtype=np.complex)
qvec=np.zeros((nqx,nqy,nqz,3),dtype=np.float)
for i in xrange(nqx):
  print i
  for j in xrange(nqy):
    for k in xrange(nqz):
      qvec[i,j,k]=np.array([float(i)/(a*nqx),float(j)/(b*nqy),float(k)/(c*nqz)]) # the q grid has to be rescaled with lattice spacing
      Vq[i,j,k,:]=np.dot(np.exp(2*1j*np.pi*np.dot(r_ext,np.array([qvec[i,j,k]]).transpose())).transpose(),U_ext)

Vq_sum=np.sum(Vq,axis=(0,1,2))/(nqx*nqy*nqz)
print 'vq summation',Vq_sum
Vq=Vq-Vq_sum


print 'writing output to',outfilename_vq
outfile_vq=h5py.File(outfilename_vq,'w')
outfile_vq.create_group('.axes')
outfile_vq['.axes/Q-points']=qvec.reshape((-1,3))
for igr,gr in enumerate(groups):
  outfile_vq[gr]=Vq[...,igr].reshape(-1)
outfile.close()

if options.show_plots:
  f4=plt.figure(4)
  plt.plot(Vq[:,0,0,plotcomp].real,linestyle='None',marker='.')
#  plt.plot(Vq[:,0,0,plotcomp].imag,linestyle='None',marker='.')
  plt.xlabel('qx')
  plt.ylabel('V(qx,0,0)')
  plt.savefig('Vq.png',dpi=300)
  f4.show()

print Vq[:,0,0,0]

print 'done, press Enter to exit'


if options.show_plots:
  raw_input()
