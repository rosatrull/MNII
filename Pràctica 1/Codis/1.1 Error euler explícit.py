#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 15:59:00 2024

@author: rosatrullcastanyer
"""

#Euler explicit
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss
from mpl_toolkits import mplot3d

#EXPLICIt

# Paràmetres inicials
n = 101  #Nombre de punts d'espai
X = np.linspace(0, 1, n)  # Mallat espacial (adimensional)
dx = 1 / (n - 1)  # Amplada del mallat espacial
t_a = 0.025  #Temps final
T_0=(0.02**2*944000)/(0.56)
gammas = [0.49, 0.25]  # Valors de gamma

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

#Càlcul de l'error relatiu
#comparació

#gamma 0,49
temperatures1 = solucions_num[0][0]
temperatures1_array = np.array(temperatures1)
solucions1_array = np.array(T_exacte)  
#err_abs = (T_exacte - temperatures)
# Operació element a element
err_abs1 = np.abs((temperatures1_array - solucions1_array))
print(err_abs1)

#gamma 0,25
temperatures2 = solucions_num[1][0]
temperatures2_array = np.array(temperatures2)
solucions2_array = np.array(T_exacte) 
#err_abs = (T_exacte - temperatures)
# Operació element a element
err_abs2 = np.abs((temperatures2_array - solucions2_array))
print(err_abs2)

#PLOT
plt.plot(z,err_abs1,label="gamma=0.49")
plt.plot(z,err_abs2,label="gamma=0,25")
plt.xlabel("z (posició)")
plt.ylabel("Error absolut")
plt.legend()
plt.grid()
plt.show()
