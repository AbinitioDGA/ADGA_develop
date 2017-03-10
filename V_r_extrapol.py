#!/usr/bin/env python
'''
This Skript reads a text file that contains r-vectors and U-values in the following format:

x   y   z   U_1   U_2   ...   U_ncomp

ncomp has to be specified as a command line option --ncomp.

Currently, also the number nd of orbitals has to be given as command line option (--nd),
but this is likely to be dropped in the next version. It is necessary only to print it 
in the output file, which is then read by Abinitio_DGA program.

Further it is necessary to provide input and output filenames,
and the number of r points for the extrapolated U(r) in x-, y- and z-direction.

The type of the fit function cannot be specified interactively, but only in the script,
approx. line 100:
                 fitfun=[yukawa | simple_exp]


sample execution line:
./V_r_extrapol.py --ncomp=6 --nd=3 ~/Data/SVO/U-00RR-vec U_r_extrapol.dat 10 10 10



Author: Josef Kaufmann,
last change: Dec 9, 2016, Ekb
'''

import numpy as np
from scipy.optimize import curve_fit

from matplotlib import pyplot as plt


import sys
import optparse



__version__ = "1.0"


parser = optparse.OptionParser(usage = "%prog [OPTIONS] input output nx ny nz",
                               version = ".".join(map(str,__version__)))

parser.add_option("--ncomp", dest="ncomp",default=0,
                 help="number of U-matrix components contained in the input file. (i.e. number of columns minus 3)")

parser.add_option("--nd", dest="nd",default=0, 
                 help="number of orbitals")


options, args = parser.parse_args()
if len(args) != 5:
  parser.error("expected input and output filename and number of r points in each direction ")

filename_in=args[0]
filename_out=args[1]
nx=int(args[2])
ny=int(args[3])
nz=int(args[4])
nr=(2*nx+1)*(2*ny+1)*(2*nz+1)
ncomp=int(options.ncomp) # number of components of U(r)
nd=int(options.nd) # number of orbitals just needed to write output file.

print nd,'Bands'
print 'U matrix has',ncomp,'components'


# read the input data
infile=open(filename_in,'r')
data=infile.readlines()
n = len(data)
r_vec=np.zeros((n,3),dtype=np.float_)
r_abs=np.zeros((n),dtype=np.float_)
U_r_in=np.zeros((n,ncomp),dtype=np.float_)
for i in xrange(n):
  r_vec[i]=data[i].split()[:3]
  r_abs[i]=np.linalg.norm(r_vec[i])
  U_r_in[i]=np.float_(data[i].split()[-ncomp:])
infile.close()

# plot the input data vs abs(r)
plotcomp=0 # component to plot
f1=plt.figure(1)
plt.plot(r_abs,U_r_in[:,plotcomp].real,marker='.',linestyle='None')
f1.show()

# do the fit
# since the function we use cannot do a constrained optimization,
# I add a "penalization" if b<=0.
# I very much hope that this is correct.
def yukawa(r,a,b):
  penalization = 0.
  if b<=0.:
    penalization=1.e20
  return a*np.exp(-b*np.abs(r))/np.abs(r) + penalization

def simple_exp(r,a,b):
  penalization = 0.
  if b<=0.:
    penalization=1.e20
  return a*np.exp(-b*np.abs(r)) + penalization
#  return 1.95*np.exp(-b*np.abs(r))


fitfun=yukawa

if fitfun==yukawa:
# if we fit to yukawa, take out the value at r=0, because we cannot fit it otherwise.
  ydata=np.zeros((n-1,ncomp),dtype=np.float_)
  rdata=np.zeros((n-1),dtype=np.float_)
  j=0
  for i in xrange(n):
    if r_abs[i] >= 0.001:
      ydata[j]=U_r_in[i]
      rdata[j]=r_abs[i]
      j=j+1
else:
  ydata=np.zeros((n,ncomp),dtype=np.float_)
  rdata=np.zeros((n),dtype=np.float_)
  j=0
  for i in xrange(n):
    ydata[j]=U_r_in[i]
    rdata[j]=r_abs[i]
    j=j+1

fitparam=np.zeros((ncomp,2),dtype=float)
fitcov = np.zeros((ncomp,2,2),dtype=float)
for iu in xrange(ncomp):
  fitparam[iu], fitcov[iu] = curve_fit(fitfun,rdata,ydata[:,iu])
  print 'component',iu
  print 'opt. fit params.',fitparam[iu]
  print 'covariance',fitcov[iu]

f2=plt.figure(2)
plt.plot(np.linspace(min(rdata),max(rdata)),fitfun(np.linspace(min(rdata),max(rdata)),fitparam[plotcomp,0],fitparam[plotcomp,1]))
plt.plot(r_abs,U_r_in[:,plotcomp].real,marker='.',linestyle='None')
f2.show()


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
  for j in xrange(2*ny+1):
    for k in xrange(2*nz+1):
      r_ext[i,j,k,:]=[(i-nx)*a,(j-ny)*b,(k-nz)*c]
      for iu in xrange(ncomp):
        U_ext[i,j,k,iu]=fitfun(np.linalg.norm(r_ext[i,j,k]),fitparam[iu,0],fitparam[iu,1])
      if np.linalg.norm(r_ext[i,j,k,:]) < np.max(r_abs):
        for l in xrange(n):
          if np.allclose(r_vec[l],r_ext[i,j,k],atol=np.min([a,b,c])/10):
            U_ext[i,j,k]=U_r_in[l]
            break
# 
      outfile.write(('{:08f}\t{:08f}\t{:08f}\n'+ncomp*'{:14f}\t'+'\n').format(r_ext[i,j,k,0],r_ext[i,j,k,1],r_ext[i,j,k,2],*U_ext[i,j,k]))
outfile.close()


f3=plt.figure(3)
plt.plot(r_ext[:,ny,nz,0],U_ext[:,ny,nz,plotcomp],marker='.',linestyle='None',label='U(x)')
# The plot in z dirction is done separately, because in HgBa2CuO4 we have different lattice spacing
# and I want to show that the fit works also in z direction. (So it is actually unnecessary for SVO)
#plt.plot(r_ext[nx,ny,:,2],U_ext[nx,ny,:,plotcomp],marker='*',linestyle='None',label='U(z)')
plt.legend()
plt.xlim(-10,10)
plt.xlabel('x,z')
f3.show()

raw_input()
