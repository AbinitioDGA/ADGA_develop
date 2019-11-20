#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.optimize
import sys

# with h5py.File('../../2dhubconv/U4_b4/output-full/adga-20191116-034405.173-output.hdf5','r') as f:
with h5py.File('../../sr2ruo4/output-full/adga-20190923-121108.189-output.hdf5','r') as f:
  siwk = f['selfenergy/nonloc/dga'][()]
  hk      = f['input/hk'][()]
  dc      = f['input/dc'][()]
  dc      = np.mean(dc,axis=1) # spins
  mudmft  = f['input/mu'][()]
  beta    = f['input/beta'][()]
  occdmft = f['input/n_dmft'][()] # required numer of electrons
  ndmft   = np.sum(occdmft) * 2
  iwfdmft = f['input/iwmax'][()] # positive number of dmft frequencies

# procedure for this whole thing
# 1) construct diagonal elements of G
# 2) extrapolate Gdiagonal towards infinite frequencies so that we can calculate n(mu)
# -> mu
# 3) extrapolate selfenergy to n4iwf + n4iwb + x frequencies (x=0)
# -> selfenergy and mu input for ADGA


def fitselfenergy(siwk, beta, fithartree, fitasymp, fext):
  '''
  Asymptotically fit the provided self-energy
  Input format is ADGA output format i.e.
  siwk.shape = ndim, ndim, kx, ky, kz, 2*iwf

  fithartree and fitasymp are the number of frequencies used
  for the asymptotic fit
  fithartree for diagonal real
  fitasymp for imaginary parts and off-diagonal real

  fext is the new number of positive freqencies

  returns siwkextended
  with siwkextended.shape = ndim, ndim, kx, ky, kz, 2*(fest)
  '''

  if fithartree < 0:
    raise ValueError('fithartree must be positive')
  if fitasymp < 0:
    raise ValueError('fitasymp must be positive')

  def fhartree(x,a,b,c):
    ''' fit function for diagonal real part '''
    return a + b/x**2 + c/x**4

  def fasymp(x,a,b,c):
    ''' fit function for off-diagonal real part and imaginary parts '''
    return a/abs(x) + b/abs(x)**2 + c/abs(x)**3

  shape         = list(siwk.shape) # shape returns a tuple
  iwf           = shape[-1]//2
  ndim          = shape[0]
  kx            = shape[2]
  ky            = shape[3]
  kz            = shape[4]
  newshape      = shape[:] # assign value not reference
  newshape[-1]  = 2*fext

  if fext < iwf:
    raise ValueError('fext must be greater than iwf')

  iwfplus = fext - iwf

  # define matsubara axes
  fmats    = np.linspace(-(iwf*2-1)*np.pi/beta,(iwf*2-1)*np.pi/beta,2*iwf) # original matsubara axis
  fmatsext = np.linspace(-((iwf+iwfplus)*2-1)*np.pi/beta,((iwf+iwfplus)*2-1)*np.pi/beta,2*(iwf+iwfplus)) # extended axis

  siwkext = np.zeros(newshape, dtype=np.complex128)
  siwkext[...,iwfplus:iwfplus+2*iwf] = siwk # center frequencies

  if fithartree > iwf:
    fithartree = iwf

  class BreakIt(Exception): pass
  # first we fit the diagonal elements
  # we have to find initial guesses so the problem is better defined
  pinit = []
  for iorb in range(ndim):
    try:
      for ikx in range(kx):
        for iky in range(ky):
          for ikz in range(kz):
            try:
              popt, _ = scipy.optimize.curve_fit(fhartree,fmats[-fithartree:],siwk[iorb,iorb,ikx,iky,ikz,-fithartree:].real)
            except RuntimeError:
              print('failed for {} {}'.format(ikx,iky))
            else:
              pinit.append(popt)
              raise BreakIt
    except BreakIt:
      continue
      print('found initial guess for hartree term for orbital {}'.format(iorb))
      print(pinit[iorb])
    else:
      raise Exception('could not find initial guess for orbital {}'.format(iorb))


  # here we only  need a couple of data points at the end to perform the fit properly
  for iorb1 in range(ndim):
    for iorb2 in range(ndim):
      print(iorb1, iorb2)
      for ikx in range(kx):
        for iky in range(ky):
          for ikz in range(kz):

            if iorb1==iorb2: # diagonal
              poptimag, _ = scipy.optimize.curve_fit(fasymp,fmats[-fitasymp:],siwk[iorb1,iorb2,ikx,iky,ikz,-fitasymp:].imag)
              poptreal, _ = scipy.optimize.curve_fit(fhartree,fmats[-fithartree:],siwk[iorb1,iorb2,ikx,iky,ikz,-fithartree:].real, p0=pinit[iorb1])

              # extend with the fitfunction and make it continous at the contact point
              siwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf:] = fhartree(fmatsext[-iwfplus:], *poptreal) \
                                                              + 1j * fasymp(fmatsext[-iwfplus:], *poptimag) \
                                                              - fhartree(fmats[-1], *poptreal) - 1j * fasymp(fmats[-1], *poptimag) \
                                                              + siwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf-1]
                                                              # last two terms to get a coninuous function
                                                              # i.e. remove any offset between the fit evaluated at the last data point of the original function
                                                              # and the original function itself

              # negative frequency extension via symmetry
              siwkext[iorb1,iorb2,ikx,iky,ikz,:iwfplus]       = np.conj(siwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf:][::-1])

            else: # offdiagonal
              poptimagpos, _ = scipy.optimize.curve_fit(fasymp,fmats[-fitasymp:],siwk[iorb1,iorb2,ikx,iky,ikz,-fitasymp:].imag)
              poptrealpos, _ = scipy.optimize.curve_fit(fasymp,fmats[-fitasymp:],siwk[iorb1,iorb2,ikx,iky,ikz,-fitasymp:].real)
              poptimagneg, _ = scipy.optimize.curve_fit(fasymp,fmats[:fitasymp],siwk[iorb1,iorb2,ikx,iky,ikz,:fitasymp].imag)
              poptrealneg, _ = scipy.optimize.curve_fit(fasymp,fmats[:fitasymp],siwk[iorb1,iorb2,ikx,iky,ikz,:fitasymp].real)

              siwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf:] = fasymp(fmatsext[-iwfplus:], *poptrealpos) \
                                                              + 1j * fasymp(fmatsext[-iwfplus:], *poptimagpos) \
                                                              - fasymp(fmats[-1], *poptrealpos) - 1j * fasymp(fmats[-1], *poptimagpos) \
                                                              + siwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf-1]

              siwkext[iorb1,iorb2,ikx,iky,ikz,:iwfplus] = fasymp(fmatsext[:iwfplus], *poptrealneg) \
                                                        + 1j * fasymp(fmatsext[:iwfplus], *poptimagneg) \
                                                        - fasymp(fmats[0], *poptrealneg) - 1j * fasymp(fmats[0], *poptimagneg) \
                                                        + siwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus]

  return siwkext



# fitselfenergy(siwk, beta=60, fithartree=30, fitasymp=15, iwfplus=140)
# sys.exit()

# def fitgreensfunction(giwk, beta, fitdiag, fitoffdiag, fext):
#   '''
#   Asymptotically fit the provided greens function
#   Input format is ADGA output format i.e.
#   giwk.shape = ndim, ndim, kx, ky, kz, 2*iwf

#   (which is constructed via the ADGA self energy)

#   fitdiag and fitoffdiag are the number of frequencies used
#   for the asymptotic fit

#   fitdiag for diagonal parts (imag and real)
#   fitoffidag for off-diagonal parts (imag and real)

#   fext is the new number of positive freqencies

#   returns siwkextended
#   with siwkextended.shape = ndim, ndim, kx, ky, kz, 2*(fext)
#   '''

#   if fitdiag < 0:
#     raise ValueError('fitdiag must be positive')
#   if fitoffdiag < 0:
#     raise ValueError('fitoffidag must be positive')

#   def fasymp(x,a,b,c):
#     ''' fit function for asymptotic behavior '''
#     return a/abs(x) + b/abs(x)**2 + c/abs(x)**3

#   shape         = list(giwk.shape) # shape returns a tuple
#   iwf           = shape[-1]//2
#   ndim          = shape[0]
#   kx            = shape[2]
#   ky            = shape[3]
#   kz            = shape[4]
#   newshape      = shape[:] # assign value not reference
#   newshape[-1]  = 2*fext

#   if fext < iwf:
#     raise ValueError('fext must be greater than iwf')

#   iwfplus = fext - iwf

#   # define matsubara axes
#   fmats    = np.linspace(-(iwf*2-1)*np.pi/beta,(iwf*2-1)*np.pi/beta,2*iwf) # original matsubara axis
#   fmatsext = np.linspace(-((iwf+iwfplus)*2-1)*np.pi/beta,((iwf+iwfplus)*2-1)*np.pi/beta,2*(iwf+iwfplus)) # extended axis

#   giwkext = np.zeros(newshape, dtype=np.complex128)
#   giwkext[...,iwfplus:iwfplus+2*iwf] = siwk # center frequencies

#   # if we provide a larger box than we have
#   # use the maximum size
#   if fitdiag > iwf:
#     fitdiag = iwf
#   if fitoffdiag > iwf:
#     fitoffdiag = iwf

#   class BreakIt(Exception): pass

#   # first we fit the diagonal elements
#   # we have to find initial guesses so the problem is better defined
#   # pinit = []
#   # for iorb in range(ndim):
#   #   try:
#   #     for ikx in range(kx):
#   #       for iky in range(ky):
#   #         for ikz in range(kz):
#   #           try:
#   #             popt, _ = scipy.optimize.curve_fit(feven,fmats[-fitdiag:],giwk[iorb,iorb,ikx,iky,ikz,-fitdiag:].imag)
#   #           except RuntimeError:
#   #             print('failed for {} {}'.format(ikx,iky))
#   #           else:
#   #             pinit.append(popt)
#   #             raise BreakIt
#   #   except BreakIt:
#   #     continue
#   #     print('found initial guess for hartree term for orbital {}'.format(iorb))
#   #     print(pinit[iorb])
#   #   else:
#   #     raise Exception('could not find initial guess for orbital {}'.format(iorb))


#   # here we only  need a couple of data points at the end to perform the fit properly
#   for iorb1 in range(ndim):
#     for iorb2 in range(ndim):
#       print(iorb1, iorb2)
#       for ikx in range(kx):
#         for iky in range(ky):
#           for ikz in range(kz):

#             if iorb1==iorb2: # diagonal
#               poptimag, _ = scipy.optimize.curve_fit(fasymp,fmats[-fitdiag:],giwk[iorb1,iorb2,ikx,iky,ikz,-fitdiag:].imag)
#               poptreal, _ = scipy.optimize.curve_fit(fasymp,fmats[-fitdiag:],giwk[iorb1,iorb2,ikx,iky,ikz,-fitdiag:].real) # , p0=pinit[iorb1])

#               # extend with the fitfunction and make it continous at the contact point
#               giwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf:] = fasymp(fmatsext[-iwfplus:], *poptreal) \
#                                                               + 1j * fasymp(fmatsext[-iwfplus:], *poptimag) \
#                                                               - fasymp(fmats[-1], *poptreal) - 1j * fasymp(fmats[-1], *poptimag) \
#                                                               + giwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf-1]
#                                                               # last two terms to get a coninuous function
#                                                               # i.e. remove any offset between the fit evaluated at the last data point of the original function
#                                                               # and the original function itself

#               # negative frequency extension via symmetry
#               giwkext[iorb1,iorb2,ikx,iky,ikz,:iwfplus]       = np.conj(giwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf:][::-1])

#             else: # offdiagonal
#               poptimagpos, _ = scipy.optimize.curve_fit(fasymp,fmats[-fitoffdiag:],giwk[iorb1,iorb2,ikx,iky,ikz,-fitoffdiag:].imag)
#               poptrealpos, _ = scipy.optimize.curve_fit(fasymp,fmats[-fitoffdiag:],giwk[iorb1,iorb2,ikx,iky,ikz,-fitoffdiag:].real)
#               poptimagneg, _ = scipy.optimize.curve_fit(fasymp,fmats[:fitoffdiag],giwk[iorb1,iorb2,ikx,iky,ikz,:fitoffdiag].imag)
#               poptrealneg, _ = scipy.optimize.curve_fit(fasymp,fmats[:fitoffdiag],giwk[iorb1,iorb2,ikx,iky,ikz,:fitoffdiag].real)

#               giwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf:] = fasymp(fmatsext[-iwfplus:], *poptrealpos) \
#                                                               + 1j * fasymp(fmatsext[-iwfplus:], *poptimagpos) \
#                                                               - fasymp(fmats[-1], *poptrealpos) - 1j * fasymp(fmats[-1], *poptimagpos) \
#                                                               + giwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus+2*iwf-1]

#               giwkext[iorb1,iorb2,ikx,iky,ikz,:iwfplus] = fasymp(fmatsext[:iwfplus], *poptrealneg) \
#                                                         + 1j * fasymp(fmatsext[:iwfplus], *poptimagneg) \
#                                                         - fasymp(fmats[0], *poptrealneg) - 1j * fasymp(fmats[0], *poptimagneg) \
#                                                         + giwkext[iorb1,iorb2,ikx,iky,ikz,iwfplus]

#   return giwkext


def construct_g(beta, mu, hk, dc, siwk):
  '''
  Construct the Greens function from its components
  G**(-1) = (iw + mu) - hk - siwk - dc

  hk, dc and siwk in the format provided from ADGA

  returns giwk
  '''

  shape = list(siwk.shape) # shape returns a tuple
  iwf   = shape[-1]//2
  ndim  = shape[0]
  fmats = np.linspace(-(iwf*2-1)*np.pi/beta,(iwf*2-1)*np.pi/beta,2*iwf)

  # ndim, ndim, kx, ky, kz, 2*iwf
  giwkinv = np.zeros(shape, dtype=np.complex128)
  giwkinv += -hk[...,None] - siwk
  giwkinv[np.arange(ndim),np.arange(ndim),...] += 1j*fmats+mu - dc[np.arange(ndim),None,None,None,None]

  giwk = np.linalg.inv(giwkinv.transpose(2,3,4,5,0,1))
  return giwk.transpose(4,5,0,1,2,3)

def calculate_occ(giwk, beta, fext, fitasymp):
  '''
  The provided Greens function G is extrapolated
  to fext positive frequencies via an asymptotic fit
  The frequency summation results in the occupation

  returns the orbital and spin resolved occupation n_{i \sigma}
  '''

  shape = list(giwk.shape) # ndim ndim 2*iwf
  ndim = shape[0]
  iwf  = shape[-1]//2

  if 2*iwf != shape[-1]:
    raise ValueError('Greens function must have an even number of matsubara frequencies')
  if fext < iwf:
    raise ValueError('fext must be larger than iwf of giw')
  if fitasymp < 0:
    raise ValueError('fitasymp must be positive')
  if fitasymp > iwf:
    fitasymp = iwf


  # positive matsubara frequencies
  fmats    = np.linspace(np.pi/beta,(iwf*2-1)*np.pi/beta,iwf)
  fmatsext = np.linspace(np.pi/beta,(fext*2-1)*np.pi/beta,fext)

  giw = np.mean(giwk.real, axis=(2,3,4)) # k mean
  # we only need the real part, the diagonal elemnts and the positive frequencies for nocc
  giwext = np.zeros((ndim,fext), dtype=np.float64)
  giwext[:,:iwf] = giw[np.arange(ndim),np.arange(ndim),iwf:]

  def fasymp(x,a,b):
    '''
    we only need to extrapolate the real part of the Greens function
    which is even
    '''
    return a/x**2 + b/x**4

  occ = []
  for iorb in range(ndim):
    popt, _ = scipy.optimize.curve_fit(fasymp, fmats[-fitasymp:], giw[iorb,iorb,-fitasymp:].real)

    giwext[iorb,iwf:] = fasymp(fmatsext[iwf:], *popt) \
                      - fasymp(fmats[-1], *popt) + giw[iorb,iorb,-1] # to make function continuous

    # factor 2 because of symmetric function ( we only sum the positive frequencies here)
    # additive 0.5 factor from the convergence factor
    occ.append(np.sum(giwext[iorb,:]) / beta * 2. + 0.5)

  occ = np.array(occ, dtype=np.float64)
  return occ



def ndeviation(mu, nrequired, beta, hk, dc, siwk, fext, fitasymp):
  giwk = construct_g(beta, mu, hk, dc, siwk)
  occ  = calculate_occ(giwk, beta, fext, fitasymp)
  print('current mu: {} - current occ: {}'.format(mu, np.sum(occ)*2))
  return nrequired - np.sum(occ)*2 # since it is spin resolved

# print(ndeviation(mu, 1.3, beta, hk, dc, siwk, 4000, 15))


ffit = 10 # number of frequencies where the fit is performed

try:
  mudga = scipy.optimize.newton(ndeviation, x0=mudmft, args=tuple((ndmft,beta,hk,dc,siwk,2*iwfdmft,ffit)), tol=1e-5)
except BaseException as e:
  print('Root finding failed')
else:
  print('Root finding succesful.\n  old mu (DMFT) = {}\n  new mu (DGA)  = {}'.format(mudmft, mudga))
  giwk = construct_g(beta, mudga, hk, dc, siwk)
  occdga  = calculate_occ(giwk, beta, 2*iwfdmft, ffit)
  print('\n  old occupation (DMFT) = {}\n  new occupation (DGA)  = {}'.format(occdmft, occdga))


# giwkext =  fitgreensfunction(giwk, beta, 15, 10, 1000)

print('Extending self-energy')
siwkext = fitselfenergy(siwk, beta, fithartree=20, fitasymp=10, fext=500)
print('Inverting extended self-energy')
giwkext = construct_g(beta, mudga, hk, dc, siwkext)

with h5py.File('test.hdf5', 'w') as h5:
  h5['siwk'] = siwkext
  h5['giwk'] = giwkext
