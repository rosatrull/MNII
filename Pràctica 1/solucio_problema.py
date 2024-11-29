import numpy as np
n = 101  # Nombre de punts d'espai
X = np.linspace(0, 1, n)  # Mallat espacial (adimensional)
dx = 1 / (n - 1)  # Amplada del mallat espacial
gamma = 0.25  # Paràmetre gamma
dt = gamma * dx**2
t_a = 0.04 #temps finals adimensional fins el que volem iterar
m = int(t_a / dt) + 1  # Nombre de passos temporals fins a arribat a ta

T_0 = (0.02**2 * 944000) / 0.56 # Valor de referència per adimensionalitzar la temperatura

L_0 = 2

C_v = 3686
P_ext = 944000
rho = 1081
t_0 = (C_v * rho * T_0) / P_ext # Valor de referència per adimensionalitzar el temps

filename=f"matriuT_expl_gamma_{gamma}.txt"
all_T = np.loadtxt(filename)
T=all_T[:,-1]
ymin=np.min(T)-273.15
ymax=np.max(T)-273.15


# Definició de la funció per verificar condicions
def verificar_condicions(T, x, t):
    z = x * L_0 #Dimensionalitzem l'espai
    regio_malalta = (z >= 0.75) & (z <= 1.25)
    regio_sana = ~regio_malalta
    
    T_C = T*T_0 - 273.15 #Passem a graus Celsius

    if np.all(T_C[regio_malalta] > 50) and np.max(T_C) < 80:
        if np.all(T_C[regio_sana] < 50):
            return True
        else:
            return False
    else:
        return False

t_optim = None


for col in range(all_T.shape[1]):  # Nombre de columnes = passos temporals +1 per incloure-ho tot
    columna = all_T[:, col]  # Seleccionar columna actual (temperatura en espai per un temps)
    t = (col * dt)*t_0  # Calcula el temps corresponent a aquesta columna

    # Verificar condicions
    if verificar_condicions(columna, X, t):
        t_optim = t  # Guardar el temps òptim si es compleixen les condicions
        break  # Sortir del bucle si trobem el primer temps que compleix

# Resultat
if t_optim is not None:
    print(f"El temps òptim és {t_optim:.2f} segons.")
else:
    print("No s'ha trobat cap temps que compleixi les condicions.")
