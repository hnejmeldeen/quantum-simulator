import numpy as np

#input potential

def V(x,y):
  if (((x*20)-10)**2+((y*20)-10)**2) == 0:
    return 0
  else:
    potential = -0.5/(((x*20)-10)**2+((y*20)-10)**2)
    if potential < 100:
      return potential
    else:
      return 100

#input wavefunction

def wavefunction(x,y):
  return 0.01
