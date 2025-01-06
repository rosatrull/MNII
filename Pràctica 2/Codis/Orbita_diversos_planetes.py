#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on M_Sn Jan  6 00:54:30 2025

@author: rosatrullcastanyer
"""

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer M_Sdern"],  # You can change the serif font here if needed
    "axes.labelsize": 14,     # Adjust as needed
    "axes.linewidth": 0.5,    # Line width for axes
    "xtick.labelsize": 12,    # Adjust tick label size
    "ytick.labelsize": 12,    
    "legend.fontsize": 10,    
    "legend.handlelength": 1.0,
    "lines.linewidth": 1,     # Width of plot lines
    "lines.markersize": 3     # Scatter dot size
})

#CONSTANTS I PARÀMETRES FÍSICS NECESSÀRIS
G=6.67428e-11                                     #cte grav universal [m^3·kg^-1·s^-1]
UA=1.496e11                                       #unitat astronòmica [m] =distància Terra-Sol
M_S=1.9885e30                                     #massa sol [kg]
T_dia_s=23.9344667*60*60                          #dia sideral [s]
T_any_s=365.25636*24*3600                         #any sideral [s]
omega=2*np.pi/T_dia_s                             #freqüència de rotació Terra [rad/s]


#NORMALITZACIÓ
r_0=UA                                              
t_0=(r_0**3/(M_S*G))**(1/2)


#DISCRETITZACIÓ
N=int(T_any_s/T_dia_s*24)                         #fem un pas cada hora
dt = T_any_s/(t_0*(N-1))                          #mida del pas temporal

#CONDICIONS INICIALS
d_af=152.1e9/r_0                                 #distància Terra-Sol a l'afeli [m]
v_af=-(29290*t_0)/r_0                            #velocitat Terra a l'afeli [m/s] !Component y!

vo = r_0 / t_0

MS = 1.9885e30 / M_S                             # massa Sol adimensional
MME = 3.3011e23/ M_S                             # massa Mercuri adimensional
MV = 4.8675e24 / M_S                             # massa Venus adimensional
MT = 5.97e24 / M_S                               # massa Terra adimensional
MM = 6.39e23 / M_S                               # massa Mart adimensional

#Valors inicials
r_Sx, r_Sy = 0, 0                                # posició inicial Sol
r_Mex, r_Mey = 69.817e9/r_0, 0                   # posició inicial Mercuri
r_Vx, r_Vy = 108.942e9/r_0, 0                    # posició inicial Venus
r_Tx, r_Ty = 152.1e9/r_0 , 0                     # posició inicial Terra
r_Mx, r_My = 249.9e9/r_0, 0                      # posició inicial Mart

v_Sx, v_Sy = 0, 0                                # velocitat inicial Sol
v_Mex, v_Mey = 0, 43599.86/vo                    # velocitat inicial Mercuri
v_Vx, v_Vy = 0, 34903.90/vo                      # velocitat inicial Venus
v_Tx, v_Ty = 0, 29290/vo                         # velocitat inicial Terra
v_Mx, v_My = 0, 22e3 /vo                         # velocitat inicial Mart

masses = [MS, MME, MV, MT, MM]
n = len(masses)  # Nombre de cossos
# llista posicions inicials
x0i = [r_Sx, r_Mex, r_Vx, r_Tx, r_Mx]
y0i = [0] * len(masses)
# llista velocitats inicials
vx0i = [0] * len(masses)
vy0i = [v_Sy, v_Mey, v_Vy, v_Ty, v_My]

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
any_M = 687 * 24 * 60 * 60                    #L'any de Mart dura 687 dies terrestres [segons]
t0 = 0
tf = any_M/ t_0
t = np.arange(t0, tf, dt)

# Llistes de posició i velocitat
xj = [[x0] for x0 in x0i]                     # Posiciones en x
yj = [[y0] for y0 in y0i]                     # Posiciones en y
vxj = [[vx0] for vx0 in vx0i]                 # Velocidades en x
vyj = [[vy0] for vy0 in vy0i]                 # Velocidades en y

# Mèt_0de de Runge-Kutta de 4t ordre (RK4)
def RK4(x, y, vx, vy, m, dt):
    xnou = []
    ynou = []
    vxnou = []
    vynou = []

    for i in range(n):

        # Càlcul de k1
        k1vx = dt * fx(x, y, m, i)
        k1vy = dt * fy(x, y, m, i)
        k1x = dt * vx[i]
        k1y = dt * vy[i]

        # Càlcul de k2
        xk = [x[j] + k1x / 2 if j == i else x[j] for j in range(n)]
        yk = [y[j] + k1y / 2 if j == i else y[j] for j in range(n)]
        k2vx = dt * fx(xk, yk, m, i)
        k2vy = dt * fy(xk, yk, m, i)
        k2x = dt * (vx[i] + k1vx / 2)
        k2y = dt * (vy[i] + k1vy / 2)

        # Càlcul de k3
        xk = [x[j] + k2x / 2 if j == i else x[j] for j in range(n)]
        yk = [y[j] + k2y / 2 if j == i else y[j] for j in range(n)]
        k3vx = dt * fx(xk, yk, m, i)
        k3vy = dt * fy(xk, yk, m, i)
        k3x = dt * (vx[i] + k2vx / 2)
        k3y = dt * (vy[i] + k2vy / 2)

        # Càlcul de k4
        xk = [x[j] + k3x if j == i else x[j] for j in range(n)]
        yk = [y[j] + k3y if j == i else y[j] for j in range(n)]
        k4vx = dt * fx(xk, yk, m, i)
        k4vy = dt * fy(xk, yk, m, i)
        k4x = dt * (vx[i] + k3vx)
        k4y = dt * (vy[i] + k3vy)

        # Nous valors de posició i velocitat
        vxnou.append(vx[i] + (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6)
        vynou.append(vy[i] + (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6)
        xnou.append(x[i] + (k1x + 2 * k2x + 2 * k3x + k4x) / 6)
        ynou.append(y[i] + (k1y + 2 * k2y + 2 * k3y + k4y) / 6)

    return xnou, ynou, vxnou, vynou

# Nous valors de posició i velocitat al llarg del temps
for i in range(len(t)-1):
    xnou, ynou, vxnou, vynou = RK4([x[-1] for x in xj], [y[-1] for y in yj],
                                                [vx[-1] for vx in vxj], [vy[-1] for vy in vyj], masses, dt)
    for i in range(n):
        xj[i].append(xnou[i])
        yj[i].append(ynou[i])
        vxj[i].append(vxnou[i])
        vyj[i].append(vynou[i])

# Gràfics de les trajectòries del planetes
fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
colors = ["yellow", "orange","green", "blue", "red"]
labels = ["Sol", "Mercuri","Venus", "Terra", "Mart"]
for i in range(n):
    ax.plot(np.array(xj[i]) * r_0, np.array(yj[i]) * r_0, label=labels[i], color=colors[i])
    ax.scatter(xj[i][0] * r_0, yj[i][0] * r_0, color=colors[i])

ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.legend(loc='upper right')

# Guardar el gràfic
plt.savefig('Orbites_N_cossos', bbox_inches='tight', dpi=300)
plt.show()
