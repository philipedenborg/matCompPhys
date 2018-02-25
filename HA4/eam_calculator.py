#!/usr/bin/env python

import numpy as np
from ase import *
from ase.calculators.eam import EAM
from scipy.interpolate import InterpolatedUnivariateSpline as spline

cutoff = 6.5

def get_calc(p):
  ''' Takes the potential parameters A, lambda, D, and 2mu and
      returns an ASE calculator object. '''
  A, lmbd, D, twomu = p

  rs   = np.linspace(0.0, cutoff, 1024)
  rhos = np.linspace(1e-10, 10, 1000)

  # Smooth cutoff function
  psi = lambda x : [(a<cutoff)*a**4/(1+a**4) for a in x]
  dpsi = lambda x : [(a<cutoff)*4*a**3/(1+a**4)**2 for a in x]

  m_phi = A*np.exp(-lmbd*rs) * psi((rs-cutoff)/3)
  m_d_phi = -lmbd*A*np.exp(-lmbd*rs) * psi((rs-cutoff)/3) + \
            A*np.exp(-lmbd*rs) * dpsi((rs-cutoff)/3) / 3
  m_embedded = -D*np.sqrt(rhos)
  m_d_embedded = -D/(2*np.sqrt(rhos))
  m_density = np.exp(-twomu*rs) * psi((rs-cutoff)/3)
  m_d_density = -twomu*np.exp(-twomu*rs) * psi((rs-cutoff)/3) + \
                 np.exp(-twomu*rs) * dpsi((rs-cutoff)/3) / 3

  m_densityf = spline(rs, m_density)
  m_d_densityf = spline(rs, m_d_density)

  m_embeddedf = spline(rhos, m_embedded)
  m_d_embeddedf = spline(rhos, m_d_embedded)

  m_phif = spline(rs, m_phi)
  m_d_phif = spline(rs, m_d_phi)

  calc = EAM(elements=['Al'],
             embedded_energy=np.array([m_embeddedf]),
             electron_density=np.array([m_densityf]),
             phi=np.array([[m_phif]]),
             d_embedded_energy=np.array([m_d_embeddedf]),
             d_electron_density=np.array([m_d_densityf]),
             d_phi=np.array([[m_d_phif]]),
             cutoff=cutoff,
             form='alloy',
             Z=[13], nr=1024, nrho=1024, dr=rs[1]-rs[0], drho=rhos[1]-rhos[0],
             lattice=['fcc'], mass=[26.982], a=[4.05])
  return calc
