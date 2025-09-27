import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import quad
from scipy.sparse import csc_array
from scipy.sparse.linalg import splu, eigsh
import importlib.util
import sys

def load_module_from_file(filepath, module_name="user_module"):
    spec = importlib.util.spec_from_file_location(module_name, filepath)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module

user_module = load_module_from_file(input("Input File? "))

n = int(float(input("Run length (seconds)? "))*30)
dt = 0.03

laser_period = 150     #time between observations  (there are 60 timesteps per second of animation)
do_check = 0

T = 0
L = int(input("Box size? "))
res = 10*L
h = 0.01
m = 0.01
count = 0
timer = 0
time_out = 1800   #how long before we assume particle got stuck behind double slits

positions = []

repeat = 1     #number of times to integrate per timestep

D = 9

X = np.linspace(0,L,res)
Y = np.linspace(0,L,res)

psi = np.zeros((res,res), dtype=complex)
for i in range(res):
  for l in range(res):
    psi[i,l] = user_module.wavefunction(X[i],Y[l])

Norm = sum(sum(np.abs(np.array(psi))**2))*(L/res)**2

psi *= 1/np.sqrt(Norm)

psi_orig = psi.copy()

V_array = np.zeros((res,res))
for i in range(res):
  for l in range(res):
    V_array[i,l] = user_module.V(i/res,l/res)


A = np.zeros((res**2,res**2),dtype=complex)

for i in range(res**2):
  A[i,i] = 1 + (dt*1j)/(2*h)*V_array[i//res,i%res] + (dt*1j*h)/(m*(L/res)**2)

for i in range(res**2 - 1):
  A[i,i+1] = -(dt*1j*h)/(4*m*(L/res)**2)
  A[i+1,i] = -(dt*1j*h)/(4*m*(L/res)**2)
  if i%res == res-1:
    A[i,i+1] = 0 + 0j
    A[i+1,i] = 0 + 0j

for i in range(res**2 - res):
  A[i,i+res] = -(dt*1j*h)/(4*m*(L/res)**2)
  A[i+res,i] = -(dt*1j*h)/(4*m*(L/res)**2)

A = csc_array(A)
lu = splu(A)

prog = []
fig = plt.figure(figsize=(10,10))
plt.xlim(0,L)
plt.ylim(0,L)
wavefunction = plt.imshow(np.abs(psi)**2,extent=(0,L,0,L),origin='lower',vmin=0,vmax=0.25)
plt.grid()
plt.colorbar()

def update(t):
  global psi, dt, count, timer, do_check, psi_orig, lu

  for k in range(repeat):
    b = []
    for i in range(res):
      for l in range(res):
        if i == 0:
          pxl = 0
        else:
          pxl = psi[i-1,l]
        if i == res-1:
          pxu = 0
        else:
          pxu = psi[i+1,l]
        if l == 0:
          pyl = 0
        else:
          pyl = psi[i,l-1]
        if l == res-1:
          pyu = 0
        else:
          pyu = psi[i,l+1]

        b.append(psi[i,l] - (dt*1j)/(2*h)*V_array[i,l]*psi[i,l] + (dt*1j*h)/(4*m*(L/res)**2)*(pxu+pyu+pxl+pyl-4*psi[i,l]))

    b = np.array(b,dtype = complex)

    psi_new = lu.solve(b)
    for i in range(res):
      for l in range(res):
        psi[i,l] = psi_new[i*res+l]

    Norm = sum(sum(np.abs(np.array(psi))**2))*(L/res)**2

    psi *= 1/np.sqrt(Norm)

    P = np.abs(psi)**2
    prob = 0
    if do_check == laser_period:
      do_check = 0
      for i in range(res):
        for l in range(res):
          if i>0.9*res:
            prob = (abs(psi[i,l])**2)*(L/res)**2
            check = np.random.rand()
            if check>prob:
              psi *= np.sqrt(1/(1-prob))
              psi[i,l] = 0 + 0j
            else:
              count += 1
              psi = psi_orig.copy()
              positions.append(l*L/res)
              timer = 0
    #do_check += 1                        #unhash this if you want laser to sample zone in the back
    timer+=1
    if timer>=time_out:
      psi = psi_orig.copy()
      timer = 0
  wavefunction.set_array(P)

  if t%(n/10)==0:
    if str(t/n*100)+'%' not in prog:
      print(str(t/n*100)+'%')
      prog.append(str(t/n*100)+'%')
  return wavefunction,

anim = FuncAnimation(fig ,update, frames=n,interval=50)
plt.rc('animation', html='html5')
anim
anim.save("animation.mp4", writer="ffmpeg", fps=30)
