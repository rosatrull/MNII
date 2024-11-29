import numpy as np
import matplotlib.pyplot as plt


T_0=(0.02**2*944000)/(0.56)

n = 101  # Nombre de punts d'espai
X = np.linspace(0, 1, n)  # Mallat espacial (adimensional)
dx = 1 / (n - 1)  # Amplada del mallat espacial

gammas = [0.25]  # Valors de gamma
L0=2 

def verificar_condicions(T, z, t):
    regio_malalta = (z >= 0.375) & (z <= 0.625) #adimensional
    regio_sana = ~regio_malalta

    if np.all(T[regio_malalta] > 50+273.15) and np.max(T) <80+273.15: #La matriu està en Kelvins per això posem les condicions en Kelvins
        if np.all(T[regio_sana] < 50+273.15):
            return True
        else:
            return False
    else:
        return False

t_optim = None

for gamma in gammas:
    filename=f"matriuT_expl_gamma_{gamma}.txt"
    all_T = np.loadtxt(filename)
    T=all_T[:,-1]
    ymin=np.min(T)-273.15
    ymax=np.max(T)-273.15
    
    dt = gamma * dx**2
    
    #mirem la columna que compleix la condicio
    for col in range(all_T.shape[1]):  # Nombre de columnes = passos temporals +1 per incloure-ho tot
        columna = all_T[:, col]  # Seleccionar columna actual (temperatura en espai per un temps)
        t = (col * dt)  # Calcula el temps corresponent a aquesta columna
        # Verificar condicions
        if verificar_condicions(columna, X, t):
            t_optim = t  # Guardar el temps òptim si es compleixen les condicions
            
            break  # Sortir del bucle si trobem el primer temps que compleix

C_v = 3686
P_ext = 944000
rho = 1081
t_0 = (C_v * rho * T_0) / P_ext # Valor de referència per adimensionalitzar el temps
#print(t_optim*t_0)         

if t_optim is not None:
    t_optim_dimensional = t_optim*t_0
    print(f"El temps òptim trobat pel mètode d'Euler implicit amb gamma = {gamma} és {t_optim_dimensional} segons")
    col_optim = int(t_optim / dt)  # Calcular la columna correponent al t_optim
    T_optim = all_T[:, col_optim]  # Obtener las temperaturas en el espacio para el t_optimo
        
    # Graficar
    fig, ax = plt.subplots(figsize=(5,4))
    ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')  
    ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
    ymin=np.min(T_optim)-273.15
    ymax=np.max(T_optim)-273.15
    #ax.plot(X*L0,np.array(T_anal)-273.15, label='Solució analítica')
    ax.scatter(X * L0, T_optim - 273.15, label=f"Temps òptim = {t_optim_dimensional:.3f}s")
    ax.hlines(50,0,2, color='black', linestyles='dashed', alpha=0.7)
    ax.set_xlabel("z (cm)")
    ax.set_ylabel(r"T "+'('+r'$\circ$'+'C)')
    ax.fill_betweenx([ymin, ymax + 2], 0.75, 1.25, color='lightcoral', alpha=0.5, edgecolor='none')
    ax.set_ylim(ymin, ymax+2)
    ax.set_xlim(0,2)
    ax.legend(loc='upper right')
    plt.savefig('Euler_expl_' + str(gamma) + '.png', bbox_inches='tight', dpi=300)
    plt.show()
