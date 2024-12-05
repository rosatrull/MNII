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
gammas = [0.51, 0.49, 0.25]  # Valors de gamma

solucions = []


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
    filename = f"matriuT_expl_gamma_{gamma}.txt"
    np.savetxt(filename, T*T_0, fmt="%.15g")



plt.figure(figsize=(10, 6))
for T_final, gamma in solucions:
    plt.plot(X, T_final, label=f"Gamma = {gamma}")
    plt.show()

'''
plt.title("Distribució de T(x, t_a) per diferents valors de gamma")
plt.xlabel("x (Posició espacial)")
plt.ylabel("T(x, t_a)")
plt.legend()
plt.grid()
plt.show()
'''


print(solucions)
