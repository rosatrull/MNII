import math
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer modern"],  # You can change the serif font here if needed
    "axes.labelsize": 14,     # Adjust as needed
    "axes.linewidth": 0.5,    # Line width for axes
    "xtick.labelsize": 12,    # Adjust tick label size
    "ytick.labelsize": 12,    
    "legend.fontsize": 10,    
    "legend.handlelength": 1.0,
    "lines.linewidth": 1,     # Width of plot lines
    "lines.markersize": 3     # Scatter dot size
})

G = 6.67428e-11  # Constant gravitacional

# Constants d'adimensionalització
Mo = 1.9885e30
xo = 1.496e11
to = np.sqrt((xo**3) / (Mo * G))
vo = xo / to

# Llista de temps
t0 = 0
tf = 365 * 24 * 60 * 60 / to
h = (60 * 60) / to
t = np.arange(t0, tf, h)

# Llei de la gravitació a l'eix X
def fx(x, y, M):
    r = np.sqrt(x*2 + y*2)
    return -(M * x) / (r**3)

# Llei de la gravitació a l'eix Y
def fy(x, y, M):
    r = np.sqrt(x*2 + y*2)
    return -(M * y) / (r**3)

# Funció de la trajectòria
def F(M,x0,y0,vx0,vy0):
  # Llistes de posició i velocitat
  xi = [x0/xo]
  yi = [y0/xo]
  vxi = [vx0/vo]
  vyi = [vy0/vo]
  # Mètode d'Euler
  for i in range(len(t) - 1):
      vxnou = vxi[i] + fx(xi[i], yi[i],M/Mo) * h
      vynou = vyi[i] + fy(xi[i], yi[i],M/Mo) * h
      xnou = xi[i] + vxnou * h
      ynou = yi[i] + vynou * h
      vxi.append(vxnou)
      vyi.append(vynou)
      xi.append(xnou)
      yi.append(ynou)

  # Redimensionalització
  xdim = [x * xo for x in xi]
  ydim = [y * xo for y in yi]
  return xdim,ydim


# Cas particular de l'òrbita el·líptica de la Terra al voltant del Sol
xel,yel= F(1.9885e30,152.101e9,0,0,29229.0)
vycirc=np.sqrt((G*1.9885e30)/(1.496e11))
xcir,ycir=F(1.9885e30,1.496e11,0,0,vycirc)

# Gràfic de la trajectòria de la Terra
fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
ax.plot(xel, yel, label="Trajectòria el·líptica", color="blue")
ax.plot(xcir, ycir, label="Trajectòria circular", color="red")
plt.scatter([0], [0], color="yellow", label="Sol")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.grid(False)
plt.legend(loc='upper right')
plt.axis("equal")
plt.show()
