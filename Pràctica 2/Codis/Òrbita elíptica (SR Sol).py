import math
import matplotlib.pyplot as plt
import numpy as np

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
    r = np.sqrt(x**2 + y**2)
    return -(M * x) / (r**3)

# Llei de la gravitació a l'eix Y
def fy(x, y, M):
    r = np.sqrt(x**2 + y**2)
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
xdim,ydim= F(1.9885e30,152.101e9,0,0,29229.0)

# Gràfic de la trajectòria de la Terra
plt.figure(figsize=(6, 6))
plt.plot(xv, yv, color="blue")
plt.scatter([0], [0], color="orange", label="Sol")  # Posición del Sol
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Òrbita de la Terra")
plt.legend()
plt.grid()
plt.axis("equal")
plt.show()
