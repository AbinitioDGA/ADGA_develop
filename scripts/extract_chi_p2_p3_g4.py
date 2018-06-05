#! /usr/bin/env python
"""
This script is meant for the comparison of chi extracted
from component-sampled G4, P3 and P2 within w2dynamics.

Mainly used to check the quality of the chi.
For that reason: Adjust the plotting commands at the
end of the script to your own likings.
"""

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from optparse import OptionParser

__version__ = 0.1
__author__ = "Matthias Pickem"

parser = OptionParser(usage = "usage: python %prog [Options]")
parser.add_option("--1pg", dest="g1", type="string",
                  help="DMFT file - one particle Green's function")
parser.add_option("--2pg", dest="g2", type="string",
                  help="DMFT file - two particle Green's function")
parser.add_option("--p3", dest="gamma", type="string",
                  help="Directly measured gamma file for comparison")
parser.add_option("--p2", dest="chi", type="string",
                  help="Directly measured gamma file for comparison")
parser.add_option("--ineq", dest="ineq", type=int, default=1,
                  help="Inequivalent atom from the DMFT calculation")
parser.add_option("--nbands", dest="nbands", type=int, default=1,
                  help="Number of correlated bands on the impurity")
parser.add_option("--bands", dest="bgroup", type=int, default=1111,
                  help="Band combination we want to consider")
parser.add_option("--spins", dest="sgroup", type=int, default=1111,
                  help="Spin combination we want to consider: gets ignored if channel options is activated")
parser.add_option("--channel", dest="channel", type=int, default=0,
                  help="Channel for the sampled group (0=spinband combination, 1=dens/magn while ignoring spins)")

def spinband2index(nbands, bands, spins):
  ''' we expect a 4 digit number for both bands and spins
      each band digit is element of [1,ndim]
      while each spin digit is lement of [1,2] = [spinup,spindown] '''
  if bands >= 1111 and bands <= (nbands*1111) and spins >= 1111 and spins <= 2222:
    b1,b2,b3,b4 = int(str(bands)[0]),int(str(bands)[1]),int(str(bands)[2]),int(str(bands)[3])
    s1,s2,s3,s4 = int(str(spins)[0]),int(str(spins)[1]),int(str(spins)[2]),int(str(spins)[3])
    index = 8*nbands**3*(2*b1+s1-3) + 4*nbands**2*(2*b2+s2-3)+2*nbands*(2*b3+s3-3)+2*b4+s4-2
    return index
  else:
    print('Wrong band- or spin-combination')
    sys.exit()

options, args = parser.parse_args()

input_g1 = h5py.File(options.g1,'r')
input_g2 = h5py.File(options.g2,'r')

beta = input_g1[".config"].attrs["general.beta"] # inverse temperature
ineq = options.ineq

giw_file = input_g1["dmft-last/ineq-{:03}/giw/value".format(ineq)][()]
# shape : nband, 2, 2*niw
niw = giw_file.shape[-1]//2
# paramagnetic solution only, average over spin component
# giw = np.zeros((nbands,2*niw), dtype=np.complex128)
giw = np.mean(giw_file, axis=1)

if (options.channel == 0):
  bandspin = spinband2index(options.nbands, options.bgroup, options.sgroup)
else:
  bandspin1 = spinband2index(options.nbands, options.bgroup, 1111)
  bandspin2 = spinband2index(options.nbands, options.bgroup, 2222)
  bandspin3 = spinband2index(options.nbands, options.bgroup, 1122)
  bandspin4 = spinband2index(options.nbands, options.bgroup, 2211)
  bandspin5 = spinband2index(options.nbands, options.bgroup, 1221)
  bandspin6 = spinband2index(options.nbands, options.bgroup, 2112)

if (options.channel == 0):
  g2_file = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin)][()]
  n4iwf = g2_file.shape[0]//2
  n4iwb = g2_file.shape[-1]//2
else:
  g2_file_1 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin1)][()]
  g2_file_2 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin2)][()]
  g2_file_3 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin3)][()]
  g2_file_4 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin4)][()]
  g2_file_5 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin5)][()]
  g2_file_6 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin6)][()]
  n4iwf = g2_file_1.shape[0]//2
  n4iwb = g2_file_1.shape[-1]//2

# this comes from user input where we start counting from 1
b1 = int(str(options.bgroup)[0])-1
b2 = int(str(options.bgroup)[1])-1
b3 = int(str(options.bgroup)[2])-1
b4 = int(str(options.bgroup)[3])-1

if (options.channel == 0):
  s1 = int(str(options.sgroup)[0])-1
  s2 = int(str(options.sgroup)[1])-1
  s3 = int(str(options.sgroup)[2])-1
  s4 = int(str(options.sgroup)[3])-1

# DATA FROM G4
# vertical
if (options.channel == 0):
  if b1 == b2 and b3 == b4 and s1 == s2 and s3 == s4:
    g2_file[:,:,n4iwb] -= beta*giw[b1,niw-n4iwf:niw+n4iwf,None]*giw[b3,None,niw-n4iwf:niw+n4iwf]

  chi_g4_sum = np.sum(g2_file, axis=(0,1))/beta**2
else:
  g2_dens = (g2_file_1 + g2_file_2 + g2_file_3 + g2_file_4)/2.0
  g2_magn = (g2_file_1 + g2_file_2 - g2_file_3 - g2_file_4 + g2_file_5 + g2_file_6 )/4.0
  if b1 == b2 and b3 == b4:
    g2_dens[:,:,n4iwb] -= 2*beta*giw[b1,niw-n4iwf:niw+n4iwf,None]*giw[b3,None,niw-n4iwf:niw+n4iwf]

  chi_g4_sum_d = np.sum(g2_dens, axis=(0,1))/beta**2
  chi_g4_sum_m = np.sum(g2_magn, axis=(0,1))/beta**2


# DATA FROM P3
input_gamma = h5py.File(options.gamma,'r')
if (options.channel == 0):
  gamma = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin)][()]
  n3iwf = gamma.shape[0]//2
  n3iwb = gamma.shape[-1]//2
else:
  gamma_1 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin1)][()]
  gamma_2 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin2)][()]
  gamma_3 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin3)][()]
  gamma_4 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin4)][()]
  gamma_5 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin5)][()]
  gamma_6 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin6)][()]
  n3iwf = gamma_1.shape[0]//2
  n3iwb = gamma_2.shape[-1]//2

occ = input_g1["dmft-last/ineq-{:03}/occ/value".format(ineq)][()]
dens = np.diagonal(np.diagonal(occ,0,-4,-2),0,-3,-2)

print(dens)
if (options.channel == 0):
  if b1 == b2 and b3 == b4 and s1 == s2 and s3 == s4:
    gamma[:,n3iwb] += beta*giw[b1,niw-n3iwf:niw+n3iwf]*(1.0-dens[b3,s3])

  chi_p3_sum = np.sum(gamma, axis=0)/beta
else:
  gamma_d = (gamma_1 + gamma_2 + gamma_3 + gamma_4)/2.0
  gamma_m = (gamma_1 + gamma_2 - gamma_3 - gamma_4 + gamma_5 + gamma_6)/4.0

  # vertical
  if b1 == b2 and b3 == b4: #and s1 == s2 and s3 == s4:
    gamma_d[:,n3iwb] += 2*beta*giw[b1,niw-n3iwf:niw+n3iwf]*(1.0-dens[b3,0])

  chi_p3_sum_d = np.sum(gamma_d, axis=0)/beta
  chi_p3_sum_m = np.sum(gamma_m, axis=0)/beta

# DATA FROM P2
input_p2 = h5py.File(options.chi,'r')
if (options.channel == 0):
  chi = input_p2["worm-001/ineq-{:03}/p2iw-worm/{:05}/value".format(ineq,bandspin)][()]
  n2iwb = chi.shape[0]//2
else:
  chi_1 = input_p2["worm-001/ineq-{:03}/p2iw-worm/{:05}/value".format(ineq,bandspin1)][()]
  chi_2 = input_p2["worm-001/ineq-{:03}/p2iw-worm/{:05}/value".format(ineq,bandspin2)][()]
  chi_3 = input_p2["worm-001/ineq-{:03}/p2iw-worm/{:05}/value".format(ineq,bandspin3)][()]
  chi_4 = input_p2["worm-001/ineq-{:03}/p2iw-worm/{:05}/value".format(ineq,bandspin4)][()]
  chi_5 = input_p2["worm-001/ineq-{:03}/p2iw-worm/{:05}/value".format(ineq,bandspin5)][()]
  chi_6 = input_p2["worm-001/ineq-{:03}/p2iw-worm/{:05}/value".format(ineq,bandspin6)][()]
  n2iwb = chi_1.shape[0]//2

if (options.channel == 0):
  if b1 == b2 and b3 == b4 and s1 == s2 and s3 == s4:
    chi[n2iwb] -= beta*(1.0-dens[b1,s1])*(1.0-dens[b3,s3])

else:
  chi_d = (chi_1 + chi_2 + chi_3 + chi_4)/2.0
  chi_m = (chi_1 + chi_2 - chi_3 - chi_4 + chi_5 + chi_6)/4.0

  # vertical
  if b1 == b2 and b3 == b4: #and s1 == s2 and s3 == s4:
    chi_d[n2iwb] -= 2*beta*(1.0-dens[b1,0])*(1.0-dens[b3,0])

if (options.channel == 0) and True:
  g = plt.figure(figsize=(12,7))
  ax1 = g.add_subplot(111)
  plt.plot(chi_g4_sum[:].real, label='G4')
  plt.plot(chi_p3_sum[n3iwb-n4iwb:n3iwb+n4iwb].real, label='P3')
  plt.plot(chi[n2iwb-n4iwb:n2iwb+n4iwb].real, label='P2')
  plt.legend(loc='best')
  plt.show()

if (options.channel == 1) and True:
  g=plt.figure(figsize=(12,10))

  ax1 = g.add_subplot(211)
  ax1.set_title(r'$\mathrm{density}$')
  plt.plot(chi_g4_sum_d[:].real, label='G4')
  plt.plot(chi_p3_sum_d[n3iwb-n4iwb:n3iwb+n4iwb].real, label='P3')
  plt.plot(chi_d[n2iwb-n4iwb:n2iwb+n4iwb].real, label='P2')
  plt.legend(loc='best')

  ax2 = g.add_subplot(212, sharex=ax1)
  ax2.set_title(r'$\mathrm{magnetic}$')
  plt.plot(chi_g4_sum_m[:].real, label='G4')
  plt.plot(chi_p3_sum_m[n3iwb-n4iwb:n3iwb+n4iwb].real, label='P3')
  plt.plot(chi_m[n2iwb-n4iwb:n2iwb+n4iwb].real, label='P2')
  plt.legend(loc='best')

  plt.show()
