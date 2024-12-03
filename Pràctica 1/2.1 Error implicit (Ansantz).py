import numpy as np
import matplotlib.pyplot as plt

#EULER EXPLICIT
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

#SOLUCIÓ ANALÍTCA    
def T_tilda(t, z, T_c, N_terms=100):
    """
    Calcula el valor de T_tilda segons la fórmula (A.7).
    
    Arguments:
    t -- temps
    z -- posició
    T_c -- temperatura de referència
    N_terms -- nombre de termes a considerar en el sumatori (precisió)
    
    Retorna:
    El valor de T_tilda(t, z)
    """
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
#err_abs = (T_exacte - temperatures)
print(len(T_exacte))
print(len(solucions[0][0]))

# Operació element a element
err_abs1 = np.abs((temperatures_array - solucions_array))
print(err_abs1)
#gamma 1
temperatures = solucions[1][0]
temperatures_array = np.array(temperatures)
solucions_array = np.array(T_exacte)  
#err_r = (T_exacte - temperatures) / temperatures
print(len(T_exacte))
print(len(solucions[1][0]))

# Operació element a element
err_abs2 = np.abs((temperatures_array - solucions_array))
print(err_abs2)

plt.plot(z,err_abs1,label="gamma=0.5")
plt.plot(z,err_abs2,label="gamma=1")
plt.xlabel("z (posició)")
plt.ylabel("Error absolut")
plt.legend()
plt.grid()
plt.show()
