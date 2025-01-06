import numpy as np
import matplotlib.pyplot as plt

# Constants
latitud = 41.9831100 * np.pi / 180  # Latitud de Girona en radians
hores_dia = 24
I0 = 10**3  

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
    return (h - 12) * 15 * np.pi / 180  

def angle_altitud(sigma, delta):
    return np.arcsin(np.cos(latitud) * np.cos(sigma) * np.cos(delta) + np.sin(latitud) * np.sin(delta))

def calcular_intensitat(dia):
    hores = np.linspace(0, hores_dia, 500)
    intensitat = []
    temps = []
    for h in hores:
        sigma = angle_horari(h)
        delta = declinacio_solar(dia, h)
        alpha_s = angle_altitud(sigma, delta)
        if alpha_s > 0:  # Només considerem quan el sol és sobre l'horitzó
            intensitat.append(I0 * np.sin(alpha_s))
            temps.append(h)
    return np.array(temps), np.array(intensitat)

def simpson(x, y):
    n = len(x) - 1
    if n % 2 != 0:  
        x = x[:-1]
        y = y[:-1]
        n = len(x) - 1
    h = (x[-1] - x[0]) / n
    integral = y[0] + y[-1] + 4 * sum(y[1:n:2]) + 2 * sum(y[2:n-1:2])
    return h * integral / 3

# Compute intensities
dia_estiu = 172  # Solstici d'estiu
dia_hivern = 355  # Solstici d'hivern

temps_estiu, intensitat_estiu = calcular_intensitat(dia_estiu)
temps_hivern, intensitat_hivern = calcular_intensitat(dia_hivern)

# Càlcul de les àrees sota les corbes amb Simpson
area_estiu_manual = simpson(temps_estiu, intensitat_estiu)
area_hivern_manual = simpson(temps_hivern, intensitat_hivern)

print(f"Àrea sota la corba (solstici d'estiu, manual Simpson): {area_estiu_manual:.2f} J/m²")
print(f"Àrea sota la corba (solstici d'hivern, manual Simpson): {area_hivern_manual:.2f} J/m²")

print(f"Energia (solstici d'estiu, manual Simpson): {area_estiu_manual * 2 * 0.2:.2f} J")
print(f"Energia (solstici d'hivern, manual Simpson): {area_hivern_manual * 2 * 0.2:.2f} J")

# Plot results
fig, ax = plt.subplots(figsize=(5,4))
ax.plot(temps_estiu, intensitat_estiu, label=f"Dia {dia_estiu} (solstici d'estiu)")
ax.plot(temps_hivern, intensitat_hivern, label=f"Dia {dia_hivern} (solstici d'hivern)")
ax.set_xlabel("Hora solar")
ax.set_ylabel("$I $"+'(W/m²)')

ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')  
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
ax.set_ylim(0,1000)

plt.grid(False)
plt.legend()
#plt.savefig('intensitat1dia.png',bbox_inches='tight', dpi=300)

plt.show()
