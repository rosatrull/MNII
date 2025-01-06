import math
import matplotlib.pyplot as plt
import numpy as np


G = 6.67428e-11  # Constant gravitacional 

# Constants d'adimensionalització
Mo = 1.9885e30   
xo = 1.496e11   
to = np.sqrt((xo**3) / (Mo * G))  
vo = xo / to                    

#Valors inicials
x0 = 1                      
y0 = 0.0                      
vx0 = 0.0                     
vy0 = np.sqrt(G * Mo / xo) / vo  

#Llista de temps 
t0 = 0
tf = 365 * 24 * 60 * 60 / to  
h = (60 * 60) / to 
t = np.arange(t0, tf, h)

#Llistes de posició i velocitat
xi = [x0]
yi = [y0]
vxi = [vx0]
vyi = [vy0]

# Llei de la gravitació
def fx(x, y, M=1):
    r = np.sqrt(x**2 + y**2)
    return -(M * x) / (r**3)

def fy(x, y, M=1):
    r = np.sqrt(x**2 + y**2)
    return -(M * y) / (r**3)

# Mètode d'Euler
for i in range(len(t) - 1):
    vxnou = vxi[i] + fx(xi[i], yi[i]) * h
    vynou = vyi[i] + fy(xi[i], yi[i]) * h
    xnou = xi[i] + vxnou * h
    ynou = yi[i] + vynou * h
    vxi.append(vxnou)
    vyi.append(vynou)
    xi.append(xnou)
    yi.append(ynou)

# Redimensionalització
xv = [x * xo for x in xi]
yv = [y * xo for y in yi]

# Gráfica de la órbita circular de la Terra al voltant del Sol
plt.figure(figsize=(6, 6))
plt.plot(xv, yv, color="blue")
plt.scatter([0], [0], color="orange", label="Sol")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Órbita circular de la Terra")
plt.legend()
plt.grid()
plt.axis("equal")
plt.show()
