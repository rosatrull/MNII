#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 19:16:09 2024

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


gammas = [0.5/2, 1/2]
n = 101  
X = np.linspace(0, 1, n)  
dx = 1 / (n - 1)  
t_a = 0.025  
T_0=(0.02**2*944000)/(0.56)

def gauss_seidel(A, b, error):
    x = np.zeros_like(b)  # Vector inicial amb zeros
    diff = 1e6   # Inicialitzar la diferència com un valor gran

    while diff > error:  # Criteri de parada
        x_old = x.copy()  # Guardar l'iteració anterior
        for i in range(len(b)):
            sum1 = np.dot(A[i, :i], x[:i])
            sum2 = np.dot(A[i, i+1:], x[i+1:])
            x[i] = (b[i] - sum1 - sum2) / A[i, i]

        diff = np.linalg.norm(x - x_old, ord=np.inf)  # Norma infinita (màxima diferència absoluta)

    return x



def matriuA(n, gamma):
    A = np.zeros((n-2, n-2)) 
    for i in range(n-2):
        A[i, i] = 1 + 2 * gamma
        if i > 0:
            A[i, i-1] = -gamma
        if i < n-3:
            A[i, i+1] = -gamma
    return A

def matriuM(n, gamma):
    M = np.zeros((n-2, n-2)) 
    for i in range(n-2):
        M[i, i] = 1 - 2 * gamma
        if i > 0:
            M[i, i-1] = gamma
        if i < n-3:
            M[i, i+1] = gamma
    return M

def solucio(gamma, n, t_a, dx):
    dt = gamma * dx**2
    m = int(t_a / dt) + 1  
    T = np.zeros((n, m)) 
    T[:, 0] = 309.65/T_0  
    T[0, :] = 309.65/T_0  
    T[-1, :] = 309.65/T_0   
    
    A = matriuA(n, gamma)
    M = matriuM(n, gamma)
    
    
    for t in range(1,m):
        b = np.dot(M, T[1:-1, t-1])  + dt
        b[0] += 2*gamma *309.65/T_0
        b[-1] += 2*gamma  *309.65/T_0
        
        T[1:-1, t] = gauss_seidel(A, b,1e-6) #Omplim les columnes interiors        
    
    return T*T_0


# Define the function that calculates T_tilda
def T_tilda(t, z, T_c, N_terms=100):
    coef = 4 / (np.pi**3)
    summation = 0
    for n in range(N_terms):
        k = 2 * n + 1
        term = (1 - np.exp(-(k**2) * (np.pi**2) * t)) / (k**3) * np.sin(k * np.pi * z)
        summation += term
    return T_c + coef * summation

# Parameters for the analytical solution
t = 0.025
z = np.linspace(0, 1, 101)
T_c = 309.65 / T_0
N_terms = 101

# Compute analytical solution
T_anal = [T_tilda(t, i, T_c) * T_0 for i in z]

L0 = 2

fig, ax = plt.subplots(figsize=(5, 4), dpi=500)
ax.plot(X * L0, np.array(T_anal) - 273.15, label='Solució analítica', color='black', linewidth=1.5)

solucions = []
for gamma in gammas:
    T = solucio(gamma, n, t_a, dx)
    solucions.append((T[:, -1], gamma))
    
ax.scatter(X * L0, solucions[0][0] - 273.15, label=f"$\\gamma = {solucions[0][1]}$", s=10)
ax.scatter(X * L0, solucions[1][0] - 273.15, label=f"$\\gamma = {solucions[1][1]}$", s=5)
# Plot horizontal line at T = 50°C
ax.hlines(50, 0, 2, color='black', linestyles='dashed', alpha=0.7)

# Shaded region between z = 0.75 cm and z = 1.25 cm
ax.fill_betweenx([np.min(T) - 273.15, np.max(T) - 273.15+2.5], 0.75, 1.25, color='lightcoral', alpha=0.5, edgecolor='none')

# Labels, legend, and axis limits
ax.set_xlabel("z (cm)")
ax.set_ylabel(r"T "+'('+r'$\circ$'+'C)')
ax.legend(loc='upper right')
ax.set_xlim(0, 2)
ax.set_ylim(np.min(T) - 273.15, np.max(T) - 273.15 + 2.5)

# Save and display the plot
plt.savefig('Crank-Nicholson', bbox_inches='tight', dpi=300)
plt.show()
"""
solucions = []
for gamma in gammas:
    T = solucio(gamma, n, t_a, dx)
    solucions.append((T[:, -1], gamma))

#print(solucions)
plt.plot(X,solucions[0][0], label=f"Gamma={gammas[0]}")
plt.plot(X,solucions[1][0], label=f"Gamma={gammas[1]}")    
plt.legend()
plt.show()
"""