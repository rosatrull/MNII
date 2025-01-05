import numpy as np

# -- DADES ---------------------------------------------------------------------------------------------
A = 2          # Àrea del panell solar (m^2)
P_max = 400    # Potència màxima (W)
I_max = 10e3   # Irradiància quan produeix potència màxima (incidència perpendicular) (W/m^2)

data = np.loadtxt("posicions_sol_placa.txt", skiprows=1) #saltem la primera línea que conté els titols
temps = data[:, 0]  # Tiempo en segundos
pos_sol = data[:, 1:]  # Posición del Sol respecto a la placa (x, y, z)

# Normal al panell (direcció z)
n = np.array([0, 0, 1])
"""
# Càlcul de theta (angle entre la incidència dels raigs solars i el vector superfície)
x = data[:, 1]
y = data[:, 2]
z = data[:, 3]

theta = np.arctan2(z, np.sqrt(x**2 + y**2)) # (radians)
"""
# -- CÀLCUL DE L'ENERGIA GENERADA -----------------------------------------------------------------------
# I = I_max * cos(theta)
energia = 0
#P_anual = 0
for i in range(len(temps) - 1):
    # Vector Sol en l'instant actual
    S = pos_sol[i]

    # Càlcul del cos de theta cos_(a^b) = a*b / |a||b|
    cos_theta = np.dot(n, S) / (np.linalg.norm(n) * np.linalg.norm(S))
    cos_theta = max(cos_theta, 0)  # Eliminar la posibilitat de tenir el sol sota l'horitzor (cos<0)

    # Potencia generada en aquest instant
    I = I_max * cos_theta
    P = P_max * (I / I_max)
    #P_anual += P
    # Pas temporal per calcular E = P * t
    delta_t = temps[i+1] - temps[i]

    # Energia generada en aquest interval
    energia += P * delta_t


print(f"Energia anual generada: {energia / 1000:.2f} kWh")
