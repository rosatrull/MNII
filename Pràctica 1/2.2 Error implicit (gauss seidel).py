import numpy as np
import matplotlib.pyplot as plt

gammas = [0.5, 1]
n = 101  
X = np.linspace(0, 1, n)  
dx = 1 / (n - 1)  
t_a = 0.025  
T_0=(0.02**2*944000)/(0.56)

def gauss_seidel(A, b): #aqui nomes estem omplint una columna, faltarà iterar sobre el temps. 
    x = np.zeros_like(b)  

    for _ in range(len(b)):  
        for i in range(len(b)):  
            sum1 = np.dot(A[i, :i], x[:i]) 
            sum2 = np.dot(A[i, i+1:], x[i+1:])  
            x[i] = (b[i] - sum1 - sum2) / A[i, i]  
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
        
        
        T[1:-1, t] = gauss_seidel(A, b) #omplim la resta de la matriu amb gauss-seidel
    
    return T*T_0


solucions = []
for gamma in gammas:
    T = solucio(gamma, n, t_a, dx)
    solucions.append((T[:, -1], gamma))  



# Definim la funció que calcula T_tilda
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

# Solució a t_a
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

print(len(T_exacte))
print(len(solucions[0][0]))

# Operació element a element
err_r1 = np.abs((temperatures_array - solucions_array) / solucions_array)
print(err_r1)

#gamma 1
temperatures = solucions[1][0]
temperatures_array = np.array(temperatures)
solucions_array = np.array(T_exacte)  
#err_r = (T_exacte - temperatures) / temperatures
print(len(T_exacte))
print(len(solucions[1][0]))

# Operació element a element
err_r2 = np.abs((temperatures_array - solucions_array) / solucions_array)
print(err_r2)

plt.plot(z,err_r1,label="gamma=0.5")
plt.plot(z,err_r2,label="gamma=1")