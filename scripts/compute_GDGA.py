#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import argparse
import sys

import numpy as np
import h5py
import scipy.optimize


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
      print('Orbital1, Orbital2: ', iorb1+1, iorb2+1)
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
    # return a/x**4 + b/x**6
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
def parse_args(args=None):
  parser = argparse.ArgumentParser(
    description='''Argument parser computation of DGA Greens function''',
    formatter_class=argparse.RawTextHelpFormatter,
    epilog="That's the end of the help")
  parser.add_argument('file',    help='ADGA output file')
  parser.add_argument('-o', '--output', help='Outputname of the HDF5 file (default="gdga.hdf5")', default='gdga.hdf5')
  parser.add_argument('--mudmft', help='Use DMFT mu instead of searching for DGA mu', default=False, action='store_true')
  parser.add_argument('--debug', help=argparse.SUPPRESS, default=False, action='store_true')

  return parser.parse_args(args)


if __name__ == '__main__':

  args = parse_args()


  try:
    with h5py.File(args.file,'r') as f:
      siwk    = f['selfenergy/nonloc/dga'][()]
      iwfdga  = siwk.shape[-1] // 2
      siw     = f['input/siw'][()]
      hk      = f['input/hk'][()]
      ndim, _, kx, ky, kz = hk.shape[:]
      dc      = f['input/dc'][()]
      dc      = np.mean(dc,axis=1) # spins
      mudmft  = f['input/mu'][()]
      beta    = f['input/beta'][()]
      occdmft = f['input/n_dmft'][()] # required numer of electrons
      ndmft   = np.sum(occdmft) * 2
      iwfdmft = f['input/iwmax'][()] # positive number of dmft frequencies
      iwbmax  = f['input/iwbmax'][()]
      iwfmax  = f['input/iwfmax'][()]

      smalliwf = iwbmax+iwfmax+20
      largeiwf = iwfdmft
  except:
    sys.exit('File is not valid.\nExiting...')

  # testing purposes: siwk from siw
  if args.debug:
    siwk = np.zeros((ndim,ndim,kx,ky,kz,2*iwfdmft), dtype=np.complex128)
    siwk[np.arange(ndim),np.arange(ndim),...] = siw[:,None,None,None,:]
    giwk = construct_g(beta, mudmft, hk, dc, siwk)
    giwk.resize(ndim,ndim,kx*ky*kz,2*iwfdmft) # this combines the k-points into the order we want
    giwk = giwk.transpose(3,2,1,0) # so we have fast access to the orbitals first
    with h5py.File('test_dmft_trunc.hdf5', 'w') as h5:
      h5['giwkext'] = giwk[iwfdmft-largeiwf:iwfdmft+largeiwf,...] # master - bubble # test
      h5['giwk']    = giwk[iwfdmft-smalliwf:iwfdmft+smalliwf,...] # everyone
    print('Done.')
    sys.exit()

  # procedure for this whole thing
  # 1) construct diagonal elements of G
  # 2) extrapolate Gdiagonal towards infinite frequencies so that we can calculate n(mu)
  # -> mu
  # 3) extrapolate selfenergy to n4iwf + n4iwb + x frequencies (x=0)
  # -> selfenergy and mu input for ADGA

  ffit =  iwfdga // 3 # number of frequencies where the fit (for the occupation) is performed)
  if args.mudmft:
    print('Enforcing mudga = mudmft')
    mudga = mudmft
    print('\n  old mu (DMFT) = {}\n  new mu (DGA)  = {}'.format(mudmft, mudga))
    giwk = construct_g(beta, mudga, hk, dc, siwk)
    occdga  = calculate_occ(giwk, beta, 2*iwfdmft, ffit)
    print('\n  old occupation (DMFT) = {}\n  new occupation (DGA)  = {}'.format(occdmft, occdga))
  else:
    print('Finding new chemical potential mudga')
    try:
      mudga = scipy.optimize.newton(ndeviation, x0=mudmft, args=tuple((ndmft,beta,hk,dc,siwk,2*iwfdmft,ffit)), tol=1e-5)
    except BaseException as e:
      print('Root finding failed')
    else:
      print('Root finding succesful.\n  old mu (DMFT) = {}\n  new mu (DGA)  = {}'.format(mudmft, mudga))
      giwk = construct_g(beta, mudga, hk, dc, siwk)
      occdga  = calculate_occ(giwk, beta, 2*iwfdmft, ffit)
      print('\n  old occupation (DMFT) = {}\n  new occupation (DGA)  = {}'.format(occdmft, occdga))

  print('Extending self-energy.')
  siwkext = fitselfenergy(siwk, beta, fithartree=iwfdga//2, fitasymp=iwfdga//4, fext=largeiwf)
  print('Constructing Greensfunction with extended self-energy.')
  giwkext = construct_g(beta, mudga, hk, dc, siwkext)
  giwkext.resize(ndim,ndim,kx*ky*kz,2*largeiwf) # this combines the k-points into the order we want
  giwkext = giwkext.transpose(3,2,1,0) # so we have fast access to the orbitals first

  print('Creating HDF5 output: {}'.format(args.output))
  with h5py.File(args.output, 'w') as h5:
    h5['giwkext'] = giwkext # master - bubble
    h5['giwk']    = giwkext[largeiwf-smalliwf:largeiwf+smalliwf,...] # everyone
  print('Done.')
