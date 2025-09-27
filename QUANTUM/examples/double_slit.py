import numpy as np

#input potential

def V(x,y):
  if 0.45<x and x<0.5:
    if y<0.37 or y>0.63 or (y>0.43 and y<0.57):
      return 0.5
    else:
      return 0
  else:
    return 0

#input wavefunction

def wavefunction(x,y):
  xinit = 2
  yinit = 5
  sigx = 10/12
  sigy = 10/8
  k0 = 8
  return np.exp(-(x-xinit)**2/(2*sigx**2)-(y-yinit)**2/(2*sigy**2))*np.cos(k0*x) + (np.exp(-(x-xinit)**2/(2*sigx**2)-(y-yinit)**2/(2*sigy**2))*np.sin(k0*x))*1j
