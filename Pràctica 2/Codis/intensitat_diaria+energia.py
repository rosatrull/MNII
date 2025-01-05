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
        if alpha_s > 0:  #només considerem quan el sol es sobre l'horitzó
            intensitat.append(I0 * np.sin(alpha_s))
            temps.append(h)
    return np.array(temps), np.array(intensitat)


dia_estiu = 172  # Solstici d'estiu
dia_hivern = 355  # Solstici d'hivern

temps_estiu, intensitat_estiu = calcular_intensitat(dia_estiu)
temps_hivern, intensitat_hivern = calcular_intensitat(dia_hivern)

# Càlcul de les àrees sota les corbes amb Simpson
area_estiu_simps = simps(intensitat_estiu, temps_estiu)
area_hivern_simps = simps(intensitat_hivern, temps_hivern)


print(f"Àrea sota la corba (solstici d'estiu, Simpson): {area_estiu_simps:.2f} J/m²")
print(f"Àrea sota la corba (solstici d'hivern, Simpson): {area_hivern_simps:.2f} J/m²")

print(f"Energia (solstici d'estiu, Simpson): {area_estiu_simps*2*0.9:.2f} J")
print(f"Energia (solstici d'hivern, Simpson): {area_hivern_simps*2*0.9:.2f} J")



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
plt.savefig('intensitat1dia.png',bbox_inches='tight', dpi=300)

plt.show()

