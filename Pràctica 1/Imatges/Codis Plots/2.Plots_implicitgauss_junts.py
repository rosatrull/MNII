#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 17:42:17 2024

@author: rosatrullcastanyer
"""

import numpy as np
import matplotlib.pyplot as plt

# Datos iniciales
gammas = [0.5, 1]
n = 101
X = np.linspace(0, 1, n)
dx = 1 / (n - 1)
t_a = 0.025
T_0 = (0.02**2 * 944000) / (0.56)
L0 = 2

# Definim la funció que calcula T_tilda
def T_tilda(t, z, T_c, N_terms=100):
    coef = 4 / (np.pi**3)
    summation = 0
    for n in range(N_terms):
        k = 2 * n + 1
        term = (1 - np.exp(-(k**2) * (np.pi**2) * t)) / (k**3) * np.sin(k * np.pi * z)
        summation += term
    return T_c + coef * summation

# Solució analítica
t = 0.025
z = np.linspace(0, 1, 101)
T_c = 309.65 / T_0
N_terms = 1000
T_exacte = [T_tilda(t, i, T_c) * T_0 for i in z]

# Gráfico con ambas gammas
fig, ax = plt.subplots(figsize=(5, 4), dpi=300)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
ax.plot(X * L0, np.array(T_exacte) - 273.15, label='Solució analítica', color='blue')

# Iterar sobre las gammas y agregar cada curva
for gamma in gammas:
    filename = f"matriuT_implGS_gamma_{gamma}.txt"
    all_T = np.loadtxt(filename)
    T = all_T[:, -1]
    dt = gamma * dx**2
    ymin = np.min(T) - 273.15
    ymax = np.max(T) - 273.15
    ax.scatter(X * L0, T - 273.15, s=10, label=f"$\\gamma = {gamma}$", alpha=0.7)

# Configuraciones del gráfico
ax.hlines(50, 0, 2, color='black', linestyles='dashed', alpha=0.7)
ax.fill_betweenx([ymin, ymax + 2], 0.75, 1.25, color='lightcoral', alpha=0.5, edgecolor='none')
ax.set_ylim(ymin, ymax + 2)
ax.set_xlim(0, 2)
ax.set_xlabel("z (cm)")
ax.set_ylabel(r"T " + '(' + r'$\circ$' + 'C)')
ax.legend(loc='upper right')

# Guardar el gráfico
plt.savefig('Euler_implGS_combined.png', bbox_inches='tight', dpi=300)
plt.show()
