#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 18:17:28 2024

@author: rosatrullcastanyer
"""

import numpy as np
import matplotlib.pyplot as plt

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


solucions = []
for gamma in gammas:
    T = solucio(gamma, n, t_a, dx)
    solucions.append((T[:, -1], gamma))

#print(solucions)
plt.plot(X,solucions[0][0]-273.15, label=f"Gamma={gammas[0]}")
plt.plot(X,solucions[1][0]-273.15, label=f"Gamma={gammas[1]}")    
plt.legend()
plt.show()