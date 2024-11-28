#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 22:27:47 2024

@author: rosatrullcastanyer
"""

import numpy as np
import matplotlib.pyplot as plt

# MÈTODE D'EUER EXPLÍCIT PER GAMMA 0.25

# Paràmetres inicials
n = 101  # Nombre de punts d'espai
X = np.linspace(0, 1, n)  # Mallat espacial (adimensional)
dx = 1 / (n - 1)  # Amplada del mallat espacial
T_0 = (0.02**2 * 944000) / 0.56
gammas = [0.51, 0.49, 0.25]  # Valors de gamma
t_a = 0.025  # Temps final (ajustarem més tard per trobar el temps òptim)

#t_0=C_v rho T_0 / P_ext
C_v = 3686
P_ext = 944000
rho = 1081

t_0 = (C_v * rho * T_0) / P_ext


# Funció per verificar les condicions
def verificar_condicions(T, z, t):
    # Condicions:
    # - Temperatura entre z=0.75 i z=1.25 > 50°C
    # - Fora d'aquest interval, per sota de 50°C
    # - Cap punt supera els 80°C
    regio_malalta = (z >= 0.75) & (z <= 1.25)
    regio_sana = ~regio_malalta

    T_C = T - 273.15  # Convertir a Celsius

    # Comprovar que tots els punts en la regió malalta estiguin per sobre de 50°C i per sota de 80°C 
    if np.all(T_C[regio_malalta] > 50) and np.max(T_C) < 80:
        # Comprovar que fora de la regió malalta tots els punts estiguin per sota de 50°C
        if np.all(T_C[regio_sana] < 50):
            return True
        else:
            return False
    else:
        return False

# Iterar sobre el temps per trobar el valor de t que compleix les condicions
t_optimo = None
#t adimensional = t / t_0
#Si volem buscar en un rang entre 30 i 90 segons que seria un resultat prou raonable -> (0.010540671179747056, 0.03162201353924116)
"""
t_i = 30/t_0
print(t_i)
t_f = 90 / t_0
print(t_f)
"""
t_values = np.linspace(0.01, 0.04, 100)  # Busquem en un interval de temps raonable

for gamma in gammas:
    dt = gamma * dx**2
    m = int(t_a / dt) + 1  # Nombre de passos temporals fins a arribat a ta
    T = np.zeros((n, m))  # Matriu T(x, t)

    # Condicions inicials
    T[:, 0] = 309.65 / T_0

    # Condicions de contorn
    T[0, :] = 309.65 / T_0
    T[-1, :] = 309.65 / T_0

    # Resolució amb Euler explícit
    for t in range(m - 1):
        for x in range(1, n - 1):
            T[x, t + 1] = T[x, t] + gamma * (T[x + 1, t] - 2 * T[x, t] + T[x - 1, t] + dx**2)
        
        # Convertir les temperatures de Kelvin a Celsius
        T_Celsius = T[:, t + 1] - 273.15
            
        # Comprovem les condicions per a cada temps
        if verificar_condicions(T[:, t], X, t):
            t_optimo = t * t_0  # Guardem el temps òptim
            break

    if t_optimo is not None:
        break

# Imprimir resultat
if t_optimo is not None:
    print(f"El valor òptim de t és: {t_optimo:.5f} segons")
else:
    print("No es va trobar un valor de t que compleixi les condicions.")

# Si es troba el temps òptim, mostrar el gràfic per a aquest temps
if t_optimo is not None:
    T_final = T[:, int(t_optimo / dt)] * T_0 - 273.15  # Temperatura final en Celsius

    plt.figure(figsize=(10, 6))
    plt.plot(X, T_final, label=f"Gamma = {gamma} a t = {t_optimo:.5f}s")
    plt.axhline(50, color='r', linestyle='--', label="50°C")
    plt.axhline(80, color='g', linestyle='--', label="80°C")
    plt.xlabel("Posició (z) [cm]")
    plt.ylabel("Temperatura (°C)")
    plt.legend()
    plt.title(f"Temperatura per a t òptim = {t_optimo:.5f}s")
    plt.grid(True)
    plt.show()
