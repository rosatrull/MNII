#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 15:59:00 2024

@author: rosatrullcastanyer
"""

#Euler explicit
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

#EXPLICIt

# Paràmetres inicials
n = 101  #Nombre de punts d'espai
X = np.linspace(0, 1, n)  # Mallat espacial (adimensional)
dx = 1 / (n - 1)  # Amplada del mallat espacial
t_a = 0.025  #Temps final
T_0=(0.02**2*944000)/(0.56)
gammas = [0.49,0.25]  # Valors de gamma

solucions = []
#t_final = []

for gamma in gammas:
    dt = gamma * dx**2
    m = int(t_a / dt) + 1  # Nombre de passos temporals fins a arribat a ta
    T = np.zeros((n, m))  # Matriu T(x, t). #n=nombre de files. m=nombre de columnes. Les files dons x1 a tots els temps i les columnes son tots els x a un cert temps. Aixo tambe es podria haver fet amb [[0 for _ in range(n)] for _ in range of m]

    # Condicions inicials
    T[:, 0] = 309.65/T_0

    # Condicions de contorn
    T[0, :] = 309.65/T_0
    T[-1, :] = 309.65/T_0


    #m-1 perque ja esta determinada la ultima posició
    for t in range(m - 1):
        for x in range(1, n - 1):
            T[x, t + 1] = T[x, t] + gamma * (T[x + 1, t] - 2 * T[x, t] + T[x - 1, t]+dx**2)


    solucions.append((T[:, -1]*T_0, gamma))  # Seleccionem la última columna
    #comprovacio = dt*m
    #t_final.append((comprovacio,gamma))
#print(t_final)

"""
#EXPLICIt

# Paràmetres inicials
n = 101  #Nombre de punts d'espai
X = np.linspace(0, 1, n)  # Mallat espacial (adimensional)
dx = 1 / (n - 1)  # Amplada del mallat espacial
t_a = 0.025  #Temps final
T_0=(0.02**2*944000)/(0.56)
gammas = [0.49, 0.25]  # Valors de gamma 0,51 no el posem perquè no dona solució i no ens interessa calcular-ne l'error

solucions_num = []


for gamma in gammas:
    dt = gamma * dx**2
    m = int(t_a / dt) + 1  # Nombre de passos temporals fins a arribat a ta
    T = np.zeros((n, m))  # Matriu T(x, t). #n=nombre de files. m=nombre de columnes. Les files dons x1 a tots els temps i les columnes son tots els x a un cert temps. Aixo tambe es podria haver fet amb [[0 for _ in range(n)] for _ in range of m]

    # Condicions inicials
    T[:, 0] = 309.65/T_0

    # Condicions de contorn
    T[0, :] = 309.65/T_0
    T[-1, :] = 309.65/T_0


    #m-1 perque ja esta determinada la ultima posició
    for t in range(m - 1):
        for x in range(1, n - 1):
            T[x, t + 1] = T[x, t] + gamma * (T[x + 1, t] - 2 * T[x, t] + T[x - 1, t]+dx**2)


    solucions_num.append((T[:, -1]*T_0, gamma))  # Seleccionem la última columna

"""
#Solució analítica

# Definim la funció que calcula T_tilda
def T_tilda(t, z, T_c, N_terms=100):
    # Constant inicial
    coef = 4 / (np.pi**3)
    # Sumatori
    summation = 0
    for n in range(N_terms):
        k = 2 * n + 1  # Índex imparell
        term = (1 - np.exp(-(k**2) * (np.pi**2) * t)) / (k**3) * np.sin(k * np.pi * z)
        summation += term  
    # Valor final
    return T_c + coef * summation

# Solució a t_a
t = 0.025  # temps fixat
z = np.linspace(0, 1, 101)  # Discretitzem l'interval [0,1]
T_c = 309.65/T_0  # Temperatura de referència
N_terms = 101  # Precisió del sumatori

# Calculem T_tilda per tots els valors de z

T_exacte=[]
for i in z:
    T_values=T_tilda(t, i, T_c)*T_0
    T_exacte.append(T_values)

"""
T_exacte1=[]
for i in z:
    T_values=T_tilda(0.025039,i,T_c)*T_0
    T_exacte1.append(T_values)
    

T_exacte2=[]
for i in z:
    T_values=T_tilda(0.025039,i,T_c)*T_0
    T_exacte2.append(T_values)
"""

#Càlcul de l'error absolut
#comparació

#gamma 0,49
temperatures1 = solucions[0][0]
temperatures1_array = np.array(temperatures1)
solucions1_array = np.array(T_exacte)  
#err_abs = (T_exacte - temperatures)
# Operació element a element
err_abs1 = np.abs((temperatures1_array - solucions1_array))


#gamma 0,25
temperatures2 = solucions[1][0]
temperatures2_array = np.array(temperatures2)
solucions2_array = np.array(T_exacte) 
#err_abs = (T_exacte - temperatures)
# Operació element a element
err_abs2 = np.abs((temperatures2_array - solucions2_array))

ymin = min(np.min(err_abs1), np.min(err_abs2))
ymax = max(np.max(err_abs1), np.max(err_abs2))

fig, ax = plt.subplots(figsize=(5, 4), dpi=500)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')

L0=2
#PLOT
ax.plot(z*L0,err_abs1,label=r"$\gamma=0.49$")
ax.plot(z*L0,err_abs2,label=r"$\gamma=0,25$")
ax.set_xlim(0, 2)
ax.set_ylim(min(np.min(err_abs1), np.min(err_abs2)),max(np.max(err_abs1), np.max(err_abs2))+0.0005)
ax.fill_betweenx([min(np.min(err_abs1), np.min(err_abs2)),max(np.max(err_abs1), np.max(err_abs2))+0.0005], 0.75 , 1.25 , color='lightcoral', alpha=0.5, edgecolor='none')
#ax.fill_betweenx([ymin, ymax + 0.0001], 0.75, 1.25, color='lightcoral', alpha=0.5, edgecolor='none')
ax.set_xlabel("z (cm)")
ax.set_ylabel("Error absolut")
ax.legend(loc="upper right")
plt.savefig('Error explícit', bbox_inches='tight', dpi=300)
plt.show()
