import numpy as np
import matplotlib.pyplot as plt

gammas = [0.5, 1]
n = 101  
X = np.linspace(0, 1, n)  
dx = 1 / (n - 1)  
t_a = 0.025  
T_0=(0.02**2*944000)/(0.56)

def gauss_seidel(A, b, error):
    x = np.zeros_like(b)  
    diff = 1e6  
    iterations = 0

    while diff > error:  # Criteri convergencia
        x_old = x.copy()  # Guardar l'iteració anterior
        for i in range(len(b)):
            sum1 = np.dot(A[i, :i], x[:i])
            sum2 = np.dot(A[i, i+1:], x[i+1:])
            x[i] = (b[i] - sum1 - sum2) / A[i, i]
        diff = np.linalg.norm(x - x_old, ord=np.inf)  # Norma infinita (màxima diferència absoluta)

    return x


def matriu(n, gamma):
    A = np.zeros((n-2, n-2)) #la matriu és en realitat una matriu nxm, però com tenim 2 condicions de contorm, tenim n-2 files. i serien m-1 columnes però em surt error si no es quadrada. Imagino q deu ser igual.
    for i in range(n-2):
        A[i, i] = 1 + 2 * gamma
        if i > 0:
            A[i, i-1] = -gamma
        if i < n-3:
            A[i, i+1] = -gamma
    return A


def solucio(gamma, n, t_a, dx):
    dt = gamma * dx**2
    m = int(t_a / dt) + 1  
    T = np.zeros((n, m))  
    T[:, 0] = 309.65/T_0  
    T[0, :] = 309.65/T_0 
    T[-1, :] = 309.65/T_0  
    
    A = matriu(n, gamma)
    
    
    for t in range(1, m):
        b = T[1:-1, t-1] + dt  #b és la columna de Temperatura a t-1, sense tenir en comte les files i columnes que ja hem ompler amb les condicions incials
        b[0] += gamma * T[0, t]  
        b[-1] += gamma * T[-1, t]  
        
        
        T[1:-1, t] = gauss_seidel(A, b,1e-8) #omplim la resta de la matriu amb gauss-seidel
    
    return T*T_0


solucions = []
for gamma in gammas:
    T = solucio(gamma, n, t_a, dx)
    solucions.append((T[:, -1], gamma))  
    filename = f"matriuT_implGS_gamma_{gamma}.txt"
    np.savetxt(filename, T, fmt="%.15g")


plt.figure(figsize=(10, 6))
for T_final, gamma in solucions:
    plt.plot(X, T_final, label=f"Gamma = {gamma}")


plt.legend()
plt.show()
