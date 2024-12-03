import numpy as np
import matplotlib.pyplot as plt

T_0=(0.02**2*944000)/(0.56)

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
N_terms = 101  # Precisió del sumatori

# Calculem T_tilda per tots els valors de z
T=[]
for i in z:
    T_values=T_tilda(t, i, T_c)*T_0
    T.append(T_values)
    
#Redfinició eix z
z=np.linspace(0,2,101)
# Fem el plot
plt.plot(z, T, label=f"t = {t}")
plt.xlabel("z [cm]")
plt.ylabel("T [K]")
plt.title("Solució analítica a t_a")
plt.fill_between(z, y1=309, y2=330, where=(z >= 0.75) & (z <= 1.25), color='red', alpha=0.5, label="Zona malalta")
plt.legend()
plt.grid()
plt.show()
