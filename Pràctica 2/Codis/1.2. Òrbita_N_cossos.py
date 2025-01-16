import math
import matplotlib.pyplot as plt
import numpy as np

#Paràmetres de format dels gràfics
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

G = 6.67428e-11  # Constant de la gravitació

# Constants d'adimensionalització
Mo = 1.9885e30
xo = 1.496e11
to = np.sqrt((xo**3) / (Mo * G))
vo = xo / to

#Valors inicials
n = 6  # Nombre de cossos
mi = [1.9885e30/Mo, 3.3011e23/Mo,4.8675e24/Mo, 5.972e24/Mo, 6.4171e23 / Mo, 1.8982e27 / Mo]  
vx0i = [0, 0, 0, 0,0,0]
x0i = [0, 69.817e9 / xo,  108.942e9 / xo, 152.101e9 / xo, 249.229e9 / xo, 816.618e9 / xo]              
y0i = [0, 0, 0, 0, 0, 0]
vy0i = [0, 43599.86 / vo, 34903.90 / vo, 29229.0 / vo, 21970.0 / vo,12444.0 / vo] 
   

# Llei de la gravitació a l'eix X
def fx(x, y, m, i):
    fxi = 0
    for k in range(n):
        if k != i:
            dx = x[k] - x[i]
            dy = y[k] - y[i]
            r = np.sqrt(dx**2 + dy**2)
            fxi += m[k] * dx / r**3
    return fxi

# Llei de la gravitació a l'eix Y
def fy(x, y, m, i):
    fyi = 0
    for k in range(n):
        if k != i:
            dx = x[k] - x[i]
            dy = y[k] - y[i]
            r = np.sqrt(dx**2 + dy**2)
            fyi += m[k] * dy / r**3
    return fyi

# Llista de temps
h = 3600/to
t0 =0
tf = (11.86*365 * 24 * 3600)/to
t = np.arange(t0, tf, h)

# Llistes de posició i velocitat
xj = [[x0] for x0 in x0i]  # Posiciones en x
yj = [[y0] for y0 in y0i]  # Posiciones en y
vxj = [[vx0] for vx0 in vx0i]  # Velocidades en x
vyj = [[vy0] for vy0 in vy0i]  # Velocidades en y

# Mètode de Runge-Kutta de 4t ordre (RK4)
def RK4(x, y, vx, vy, m, h):
    xnou = []
    ynou = []
    vxnou = []
    vynou = []

    for i in range(n):

        # Càlcul de k1
        k1vx = h * fx(x, y, m, i)
        k1vy = h * fy(x, y, m, i)
        k1x = h * vx[i]
        k1y = h * vy[i]

        # Càlcul de k2
        xk = [x[j] + k1x / 2 if j == i else x[j] for j in range(n)]
        yk = [y[j] + k1y / 2 if j == i else y[j] for j in range(n)]
        k2vx = h * fx(xk, yk, m, i)
        k2vy = h * fy(xk, yk, m, i)
        k2x = h * (vx[i] + k1vx / 2)
        k2y = h * (vy[i] + k1vy / 2)

        # Càlcul de k3
        xk = [x[j] + k2x / 2 if j == i else x[j] for j in range(n)]
        yk = [y[j] + k2y / 2 if j == i else y[j] for j in range(n)]
        k3vx = h * fx(xk, yk, m, i)
        k3vy = h * fy(xk, yk, m, i)
        k3x = h * (vx[i] + k2vx / 2)
        k3y = h * (vy[i] + k2vy / 2)

        # Càlcul de k4
        xk = [x[j] + k3x if j == i else x[j] for j in range(n)]
        yk = [y[j] + k3y if j == i else y[j] for j in range(n)]
        k4vx = h * fx(xk, yk, m, i)
        k4vy = h * fy(xk, yk, m, i)
        k4x = h * (vx[i] + k3vx)
        k4y = h * (vy[i] + k3vy)

        # Nous valors de posició i velocitat
        vxnou.append(vx[i] + (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6)
        vynou.append(vy[i] + (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6)
        xnou.append(x[i] + (k1x + 2 * k2x + 2 * k3x + k4x) / 6)
        ynou.append(y[i] + (k1y + 2 * k2y + 2 * k3y + k4y) / 6)

    return xnou, ynou, vxnou, vynou

# Nous valors de posició i velocitat al llarg del temps
for i in range(len(t)-1):
    xnou, ynou, vxnou, vynou = RK4([x[-1] for x in xj], [y[-1] for y in yj],
                                                [vx[-1] for vx in vxj], [vy[-1] for vy in vyj], mi, h)
    for i in range(n):
        xj[i].append(xnou[i])
        yj[i].append(ynou[i])
        vxj[i].append(vxnou[i])
        vyj[i].append(vynou[i])

# Gràfics de les trajectòries del planetes
fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
colors = ["yellow", "red","purple", "blue",'orange','brown']
labels = ["Sol", "Mercuri","Venus", "Terra","Mart","Júpiter"]
for i in range(n):
    if i==0:
      plt.scatter(xj[i][0] * xo, yj[i][0] * xo, color=colors[i],label=labels[i])
    else:
      plt.plot(np.array(xj[i]) * xo, np.array(yj[i]) * xo, label=labels[i], color=colors[i])
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.grid(False)
plt.legend(loc='upper right')
plt.axis("equal")
plt.show()
