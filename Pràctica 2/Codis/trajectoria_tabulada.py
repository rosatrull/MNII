import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#intentem no tocar això, així tots els plots quedaràn amb mides iguals. 

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer modern"],  # You can change the serif font here if needed
    "axes.labelsize": 14,     # Adjust as needed
    "axes.linewidth": 0.5,    # Line width for axes
    "xtick.labelsize": 12,    # Adjust tick label size
    "ytick.labelsize": 12,    
    "legend.fontsize": 10,    
    "legend.handlelength": 1.0,
    "lines.linewidth": 1,     # Width of plot lines
    "lines.markersize": 3     # Scatter dot size
})



latitud = 41.9831100 * np.pi / 180  # Latitud de Girona en radians
dies_any = 365
hores_dia = 24

# Funcions per calcular els angles
def declinacio_solar(N, h):
    x = 2 * np.pi * (N - 1 + (h - 12) / 24) / 365
    return (
        0.006918
        - 0.399912 * np.cos(x)
        + 0.070257 * np.sin(x)
        - 0.006758 * np.cos(2 * x)
        + 0.000907 * np.sin(2 * x)
        - 0.002697 * np.cos(3 * x)
        + 0.001480 * np.sin(3 * x)
    )

def angle_horari(h):
    return (h - 12) * 15 * np.pi / 180  # Angle horari solar en radians

def angle_altitud(sigma, delta):
    return np.arcsin(np.cos(latitud) * np.cos(sigma) * np.cos(delta) + np.sin(latitud) * np.sin(delta))

def angle_azimut(sigma, delta, alpha_s):
    sin_gamma = np.sin(sigma)
    cos_gamma = (np.cos(sigma) * np.sin(latitud) - np.tan(delta) * np.cos(latitud))
    gamma = np.arctan2(sin_gamma, cos_gamma)
    return gamma % (2 * np.pi)  

# Generar les coordenades 3D del camí del Sol
def cami_sol(dia):
    hores = np.linspace(0, hores_dia, 500)  # Temps des de q surt el sol fins q es pon
    x, y, z = [], [], []
    for h in hores:
        sigma = angle_horari(h)
        delta = declinacio_solar(dia, h)
        alpha_s = angle_altitud(sigma, delta)

        # Filtrar els valors d'altitud negativa (nit)
        if alpha_s < 0:
            continue

        gamma_s = angle_azimut(sigma, delta, alpha_s)

        # Convertir coordenades esfèriques a cartesianes per a la visualització
        x.append(np.cos(alpha_s) * np.sin(gamma_s))
        y.append(np.cos(alpha_s) * np.cos(gamma_s))
        z.append(np.sin(alpha_s))

    # Assegurar que el camí comenci i acabi al pla XY
    if z:
        z[0] = 0
        z[-1] = 0
    return np.array(x), np.array(y), np.array(z)

# Representació del camí del Sol
fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(111, projection='3d')

# Trajectòries del Sol per diferents dies de l'any
for dia in np.arange(0, 366, 10):  # Interval de 10 dies
    x, y, z = cami_sol(dia)
    ax.plot(x, y, z, color='black')


ax.set_title("Evolució 3D del Sol des de la sortida fins a la posta")
ax.set_xlabel("Est")
ax.set_ylabel("Nord", labelpad=10)
ax.set_zlabel("Altitud", labelpad=-30) 
ax.xaxis.pane.set_visible(False)
ax.yaxis.pane.set_visible(False)
ax.zaxis.pane.set_visible(True)  
ax.xaxis.pane.set_edgecolor('grey')
ax.yaxis.pane.set_edgecolor('grey')
ax.xaxis.pane.set_facecolor((0.9, 0.9, 0.9, 1.0))  
ax.yaxis.pane.set_facecolor((0.9, 0.9, 0.9, 1.0))  

ax.set_xticks(np.arange(-1, 1+0.5, 0.5),minor=False)  
ax.set_yticks(np.arange(-0.5, 1+0.5, 0.5), minor=False) 
ax.set_zticks(np.arange(0, 1+0.5, 0.5), minor=False)

ax.grid(False)

plt.show()
