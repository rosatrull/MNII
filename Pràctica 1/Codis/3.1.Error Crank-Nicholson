#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 19:29:41 2024

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


solucions = []
for gamma in gammas:
    T = solucio(gamma, n, t_a, dx)
    solucions.append((T[:,-1], gamma))  
    
#SOLUCIÓ ANALÍTCA    
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

# Sol analítica
t = 0.025  # temps fixat
z = np.linspace(0, 1, 101)  # Discretitzem l'interval [0,1]
T_c = 309.65/T_0  # Temperatura de referència
N_terms = 1000  # Precisió del sumatori

# Calculem T_tilda per tots els valors de z
T_exacte=[]
for i in z:
    T_values=T_tilda(t, i, T_c)*T_0
    T_exacte.append(T_values)


#Comparació
#gamma 0.5
temperatures = solucions[0][0]
temperatures_array = np.array(temperatures)
solucions_array = np.array(T_exacte)  
#err_r = (T_exacte - temperatures) / temperatures
print(len(T_exacte))
print(len(solucions[0][0]))

# Operació element a element
err_r1 = np.abs(temperatures_array - solucions_array) 
print(err_r1)
#gamma 1
temperatures = solucions[1][0]
temperatures_array = np.array(temperatures)
solucions_array = np.array(T_exacte)  
#err_r = (T_exacte - temperatures) / temperatures
print(len(T_exacte))
print(len(solucions[1][0]))

# Operació element a element
err_r2 = np.abs(temperatures_array - solucions_array) 
print(err_r2)

plt.plot(z,err_r1,label="gamma=0.5")
plt.plot(z,err_r2,label="gamma=1")
plt.xlabel("z (posició)")
plt.ylabel("T")
plt.title("Solució de T a temps fix")
plt.legend()
plt.grid()
plt.show()
