#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit

from matplotlib import pyplot as plt


import sys
import optparse



__version__ = "1.0"


parser = optparse.OptionParser(usage = "%prog [OPTIONS] <hdf-file>",
                               version = ".".join(map(str,__version__)))


options, args = parser.parse_args()
if len(args) != 5:
  parser.error("expected input and output filename and number of r points in each direction ")

filename_in=args[0]
filename_out=args[1]
nx=int(args[2])
ny=int(args[3])
nz=int(args[4])
nr=(2*nx+1)*(2*ny+1)*(2*nz+1)
nd=1


# read the input data
infile=open(filename_in,'r')
data=infile.readlines()
n = len(data)
r_vec=np.zeros((n,3),dtype=np.float_)
r_abs=np.zeros((n),dtype=np.float_)
U_r_in=np.zeros((n,nd**2,nd**2),dtype=np.float_)
for i in xrange(n):
  r_vec[i]=data[i].split()[:3]
  r_abs[i]=np.linalg.norm(r_vec[i])
  U_r_in[i,0,0]=np.float_(data[i].split()[-1])
infile.close()

# plot the input data vs abs(r)
f1=plt.figure(1)
plt.plot(r_abs,U_r_in[:,0,0].real,marker='.',linestyle='None')
f1.show()

# take out the value at r=0, because we cannot fit it otherwise.
ydata=np.zeros((n-1,nd**2,nd**2),dtype=np.float_)
rdata=np.zeros((n-1),dtype=np.float_)
j=0
for i in xrange(n):
  if r_abs[i] >= 0.001:
    ydata[j]=U_r_in[i]
    rdata[j]=r_abs[i]
    j=j+1

# do the fit
def yukawa(r,a,b):
  return a*np.exp(-b*np.abs(r))/np.abs(r)


popt, pcov = curve_fit(yukawa,rdata,ydata[:,0,0])

print popt
print pcov

f2=plt.figure(2)
plt.plot(np.linspace(min(rdata),max(rdata)),yukawa(np.linspace(min(rdata),max(rdata)),popt[0],popt[1]))
plt.plot(r_abs,U_r_in[:,0,0].real,marker='.',linestyle='None')
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
U_ext=np.zeros((2*nx+1,2*ny+1,2*nz+1,nd**2,nd**2),dtype=np.float_)


outfile=open(filename_out,'w')
outfile.write('{} {} {} {} {} \n'.format(nr,nd,a,b,c))
for i in xrange(2*nx+1):
  for j in xrange(2*ny+1):
    for k in xrange(2*nz+1):
      r_ext[i,j,k,:]=[(i-nx)*a,(j-ny)*b,(k-nz)*c]
      U_ext[i,j,k]=yukawa(np.linalg.norm(r_ext[i,j,k]),popt[0],popt[1])
      if np.linalg.norm(r_ext[i,j,k,:]) < np.max(r_abs):
        for l in xrange(n):
          if np.allclose(r_vec[l],r_ext[i,j,k],atol=np.min([a,b,c])/10):
            U_ext[i,j,k]=U_r_in[l]
            break
# writing to file should be done in nd*nd blocks, just 1*1 is trivial
      outfile.write('{:08f}\t{:08f}\t{:08f}\n{:14f}\n'.format(r_ext[i,j,k,0],r_ext[i,j,k,1],r_ext[i,j,k,2],U_ext[i,j,k,0,0]))
outfile.close()

f3=plt.figure(3)
plt.plot(r_ext[:,ny,nz,0],np.log10(U_ext[:,ny,nz,0,0]),marker='.',linestyle='None')
f3.show()

raw_input()

