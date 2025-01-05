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
    return (h - 12) * 15 * np.pi / 180  # Angle horari solar en radians

def angle_altitud(sigma, delta):
    return np.arcsin(np.cos(latitud) * np.cos(sigma) * np.cos(delta) + np.sin(latitud) * np.sin(delta))


def calcular_intensitat_diaria(dia):
    hores = np.linspace(0, hores_dia, 500)  
    intensitats = []
    for h in hores:
        sigma = angle_horari(h)
        delta = declinacio_solar(dia, h)
        alpha_s = angle_altitud(sigma, delta)
        if alpha_s > 0:  
            intensitats.append(I0 * np.sin(alpha_s))
    if intensitats:
        return np.mean(intensitats)
    return 0 


dies = np.arange(1, 366) 
intensitats_mitjanes = [calcular_intensitat_diaria(dia) for dia in dies]

area= simps(intensitats_mitjanes, dies)

# Mostrem les àrees calculades
print(f"Àrea sota la corba (solstici d'estiu, Simpson): {area:.2f} J/m²")
print(f"Energia (solstici d'hivern, Simpson): {area*2*0.9:.2f} J")
fig, ax = plt.subplots(figsize=(5,4))

ax.set_xlabel("Dia de l'any")
ax.set_ylabel("$I $"+'(W/m²)')

ax.plot(dies, intensitats_mitjanes, label="Intensitat mitjana diària")

ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')  
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
plt.axvline(172, color='orange', linestyle='--', label="Solstici d'estiu (dia 172)")
plt.axvline(355, color='blue', linestyle='--', label="Solstici d'hivern (dia 355)")

ax.set_ylim(270,600)
plt.grid(False)
plt.legend(loc='upper right')
plt.savefig('intensitat mitjana.png',bbox_inches='tight', dpi=300)

plt.show()

