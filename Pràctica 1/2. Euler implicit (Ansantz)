import numpy as np
import matplotlib.pyplot as plt

# Paràmetres inicials
n = 101  
X = np.linspace(0, 1, n) 
dx = 1 / (n - 1)  
t_a = 0.025  
gammas = [0.5, 1]  
T_0=(0.02**2*944000)/(0.56)
solucions = []


for gamma in gammas:
    dt = gamma * dx**2
    m = int(t_a / dt) + 1                               #501-gamma=0.5 // 251-gamma=1.0
    #print(m)
    T = np.zeros((n, m))                                #(files=101;columnes=501/251)
    T[:, 0] = 309.65/T_0                                #C.I. [ADIMENSIONALITZADES]
    T[0, :] = 309.65/T_0                                #C.C.
    T[-1, :] = 309.65/T_0                               #C.C.
    #print(T)
    #print(T.shape[1])
    
    for i in range(1, T.shape[1]):                      #Iterem sobre les columnes
        # Definir T^{i-1}
        T_prev = T[:, i-1]                              #El 1r cop agafa una columna sencera de 309, el 2n cop nmes 309 als extrems !!!
        T_nou_i = T_prev.copy()                         #Per la comparació
        T_nou = T_prev.copy()                           #A modificar

        error = 1e6                                     #Error gran per entrar al bucle
        while error > 1e-6:                             #Continuar fins que l'error < límit
            # CALCULEM ELS NOUS ELEMENTS
            for j in range(T.shape[0]):                 #Iterar sobre les files
                if j == 0:
                    T_nou[j] = 309.65/T_0
                elif j == n - 1:
                    T_nou[j] = 309.65/T_0
                else:
                    T_nou[j] = (T_prev[j] + gamma * (T_nou_i[j + 1] + T_nou_i[j - 1]) + dt) / (1 + 2 * gamma)
    
            # Calcul de l'error màxim (norma infinita)
            error = np.linalg.norm(T_nou - T_nou_i, ord=np.inf)  
            if error > 1e-6:
                T_nou_i = T_nou.copy()  # Actualitzar per continuar iterant
    
        # Quan convergeix, desar la columna  
        T[:, i] = T_nou
    solucions.append((T[:,-1]*T_0,gamma))              #RETORNEM A UNITATS DIMENSIONALS
    
     
#print(f"Aixo es la solució= {solucions}")
#print(len(solucions[0][0]))
#print(len(solucions[1][0]))    
plt.plot(X, solucions[0][0], label=f"Gamma={solucions[0][1]}")
plt.plot(X, solucions[0][0], label=f"Gamma={solucions[1][1]}")
plt.legend()
plt.show()            
        
