#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 08:17:34 2025

@author: rosatrullcastanyer
"""

import numpy as np
import matplotlib.pyplot as plt

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

pos_i=np.zeros((N,3))                             #(x,y,z) a t_n
pos_f=np.zeros((N,3))                             #(x,y,z) a t_n+1
pos_f[0,0]=d_af
pos_f[0,1]=0

vel_i=np.zeros((N,3))                             #(vx,vy,vz) a t_n
vel_f=np.zeros((N,3))                             #(vx,vy,vz) a t_n+1
vel_f[0,0]=0
vel_f[0,1]=v_af


#RUNGE-KUTTA 4-----------------------------------------------------------------------------------------
t=0
for j in range(1, N):
    t = t + dt

    # Actualitzem k1
    k1x = dt * vel_f[j-1][0]
    k1y = dt * vel_f[j-1][1]
    k1vx = -dt * pos_f[j-1][0] / (pos_f[j-1][0]**2 + pos_f[j-1][1]**2)**(3/2)
    k1vy = -dt * pos_f[j-1][1] / (pos_f[j-1][0]**2 + pos_f[j-1][1]**2)**(3/2)

    # Actualitzem k2
    k2x = dt * (vel_f[j-1][0] + k1vx / 2)
    k2y = dt * (vel_f[j-1][1] + k1vy / 2)
    k2vx = -dt * (pos_f[j-1][0] + k1x / 2) / (((pos_f[j-1][0] + k1x / 2)**2 + (pos_f[j-1][1] + k1y / 2)**2)**(3/2))
    k2vy = -dt * (pos_f[j-1][1] + k1y / 2) / (((pos_f[j-1][0] + k1x / 2)**2 + (pos_f[j-1][1] + k1y / 2)**2)**(3/2))

    # Actualitzem k3
    k3x = dt * (vel_f[j-1][0] + k2vx / 2)
    k3y = dt * (vel_f[j-1][1] + k2vy / 2)
    k3vx = -dt * (pos_f[j-1][0] + k2x / 2) / (((pos_f[j-1][0] + k2x / 2)**2 + (pos_f[j-1][1] + k2y / 2)**2)**(3/2))
    k3vy = -dt * (pos_f[j-1][1] + k2y / 2) / (((pos_f[j-1][0] + k2x / 2)**2 + (pos_f[j-1][1] + k2y / 2)**2)**(3/2))

    # Actualitzem k4
    k4x = dt * (vel_f[j-1][0] + k3vx)
    k4y = dt * (vel_f[j-1][1] + k3vy)
    k4vx = -dt * (pos_f[j-1][0] + k3x) / (((pos_f[j-1][0] + k3x)**2 + (pos_f[j-1][1] + k3y)**2)**(3/2))
    k4vy = -dt * (pos_f[j-1][1] + k3y) / (((pos_f[j-1][0] + k3x)**2 + (pos_f[j-1][1] + k3y)**2)**(3/2))

    # Posició terra respecte el Sol
    pos_f[j][0] = pos_f[j-1][0] + (k1x + 2*k2x + 2*k3x + k4x) / 6
    pos_f[j][1] = pos_f[j-1][1] + (k1y + 2*k2y + 2*k3y + k4y) / 6
    vel_f[j][0] = vel_f[j-1][0] + (k1vx + 2*k2vx + 2*k3vx + k4vx) / 6
    vel_f[j][1] = vel_f[j-1][1] + (k1vy + 2*k2vy + 2*k3vy + k4vy) / 6 
    
    
#DIMENSIONALITZEM
xd=pos_f[:, 0]*r_0
yd=pos_f[:, 1]*r_0
#zd=pos_f[:, 2]*r_0


#MÈTODE D'EULER------------------------------------------------------------------------------------
pos_f_e=np.zeros((N,3))                             #(x,y,z) a t_n+1
pos_f_e[0,0]=d_af
pos_f_e[0,1]=0

vel_f_e=np.zeros((N,3))                             #(vx,vy,vz) a t_n+1
vel_f_e[0,0]=0
vel_f_e[0,1]=v_af

#MÈTODE D'EULER----------------------------------------------------------------------------------------------------------
def acceleracio_x(x,y):
    ax = - x / (x**2 + y**2)**(3/2)
    return ax
def acceleracio_y(x,y):
    ay = - y / (x**2 + y**2)**(3/2)
    return ay

for j in range(1,N):

    #Actualitzem velocitats
    vel_f_e[j,0] = vel_f_e[j-1,0] + acceleracio_x(pos_f_e[j-1,0], pos_f_e[j-1,1]) * dt
    vel_f_e[j,1] = vel_f_e[j-1,1] + acceleracio_y(pos_f_e[j-1,0], pos_f_e[j-1,1]) * dt
    #Actualitzem posicions
    pos_f_e[j,0] = pos_f_e[j-1,0] + vel_f_e[j-1,0] * dt 
    pos_f_e[j,1] = pos_f_e[j-1,1] + vel_f_e[j-1,1] * dt 
   

#DIMENSIONALITZEM
xd_e=pos_f_e[:, 0]*r_0
yd_e=pos_f_e[:, 1]*r_0
#zd=pos_f[:, 2]*r_0

#SOLUCIÓ ANALÍTICA---------------------------------------------------------------------------------------
def analitica(theta):
    L = 2.66e40
    M_S = 1.9885e30 
    G = 6.67428e-11
    e = 0.0167
    M_T = 5.97e24
    denominador = G * M_S * (M_T)**2 * (1 + e*np.cos(theta))
    r = L**2 / denominador
    return r

angle_theta = np.linspace(0, 2 * np.pi, N)
r_values = analitica(angle_theta)

x_an = r_values * np.cos(angle_theta) * r_0
y_an = r_values * np.sin(angle_theta) * r_0


#GRÀFIQUEM LES DIVERSES SOLUCIONS-------------------------------------------------------------------
#Euler + RK4
fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')

ax.plot(xd, yd, label='Orbita de la Terra al voltant del Sol amb RK4', color='blue')
ax.plot(xd_e, yd_e, label='Orbita de la Terra al voltant del Sol amb Euler', color='red')
#ax.plot(x_an,y_an, label='Òrbita teòrica de la Terra al voltant del Sol', color = 'black')

ax.scatter(0, 0, s=40, label="Sol", color="orange")

#ax.set_ylim(min(yd)-1e11, max(yd)+1e11)
#ax.set_xlim(min(xd)-1e11, max(xd)+1e11)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.legend(loc='upper right')

# Guardar el gràfic
plt.savefig('Orbita_Terra_Euler+RK4', bbox_inches='tight', dpi=300)
plt.show()

#Euler + RK4 + teroica
fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')

ax.plot(xd, yd, label='Orbita de la Terra al voltant del Sol amb RK4', color='blue')
ax.plot(xd_e, yd_e, label='Orbita de la Terra al voltant del Sol amb Euler', color='red')
ax.plot(x_an,y_an, label='Òrbita teòrica de la Terra al voltant del Sol', color = 'black')

ax.scatter(0, 0, s=40, label="Sol", color="orange")

#ax.set_ylim(min(yd)-1e11, max(yd)+1e11)
#ax.set_xlim(min(xd)-1e11, max(xd)+1e11)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.legend(loc='upper right')

# Guardar el gràfic
plt.savefig('Orbita_Terra_Euler+RK4+teorica', bbox_inches='tight', dpi=300)
plt.show()

#CALCULAR L'ERROR RELATIU--------------------------------------------------------------------------------------
"""
x_error_euler = pos_f_e[:,0] - x_an
y_error_euler = pos_f_e[:,1] - y_an

x_error_kutta = pos_f[:,0] - x_an
y_error_kutta = pos_f[:,1] - y_an

error_euler = np.abs(np.sqrt(x_error_euler**2 + y_error_euler**2))  
error_rk4 = np.abs(np.sqrt(x_error_kutta**2 + y_error_kutta**2)) 
"""
r_euler = np.sqrt((pos_f_e[:, 0]*r_0)**2 + (pos_f_e[:, 1]*r_0)**2)
r_kutta = np.sqrt((pos_f[:, 0]*r_0)**2 + (pos_f[:, 1]*r_0)**2)

error_euler = np.abs(r_euler - r_values)
error_rk4 = np.abs(r_kutta - r_values)

# Punts pels ticks
pi_ticks = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
pi_labels = [r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]

# Gràfic
fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')

ax.plot(angle_theta, error_euler, label="Error absolut amb el mètode d'Euler", color='red')
ax.plot(angle_theta, error_rk4, label="Error absolut amb el mètode RK4", color='blue')

# Configurar ticks de l'eix x
ax.set_xticks(pi_ticks)
ax.set_xticklabels(pi_labels)

# Etiquetes i llegenda
ax.set_xlabel(r"$\theta$ (rad)")
ax.set_ylabel("Error absolut")
ax.legend(loc='upper right')

# Guardar i mostrar el gràfic
plt.savefig('Error_absolut', bbox_inches='tight', dpi=300)
plt.show()
#Plot
"""
fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')

ax.plot(angle_theta, error_euler, label="Error absolut amb el mètode d'Euler", color='red')
ax.plot(angle_theta, error_rk4, label="Error absolut amb el mètode RK4", color='blue')

#ax.set_ylim(min(yd)-1e11, max(yd)+1e11)
#ax.set_xlim(min(xd)-1e11, max(xd)+1e11)
ax.set_xlabel(r"$\theta$ (rad))")
ax.set_ylabel("Error absolut")
ax.legend(loc='upper right')

# Guardar el gràfic
plt.savefig('Error_absolut', bbox_inches='tight', dpi=300)
plt.show()
"""
