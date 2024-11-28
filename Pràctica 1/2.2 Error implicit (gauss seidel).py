import numpy as np
import matplotlib.pyplot as plt

gammas = [0.5, 1]
n = 101  
X = np.linspace(0, 1, n)  
dx = 1 / (n - 1)  
t_a = 0.025  
T_0=(0.02**2*944000)/(0.56)
L0=2


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




for gamma in gammas:
    filename=f"matriuT_implGS_gamma_{gamma}.txt"
    all_T = np.loadtxt(filename)
    T=all_T[:,-1]
    dt = gamma * dx**2
    m = int(t_a / dt) + 1  
    ymin=np.min(T)-273.15
    ymax=np.max(T)-273.15

    #volem fer una evolució temporal, pel que a cada plot cal plotejar cada columna. 
    fig, ax = plt.subplots(figsize=(5,4))
    ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')  
    ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
    ymin=np.min(T)-273.15
    ymax=np.max(T)-273.15
    ax.plot(X*L0,np.array(T_exacte)-273.15, label='Solució analítica')
    ax.scatter(X * L0, T - 273.15, label=f"$\\gamma = {gamma}$\n$\\hat{{t}}=0.025$")
    ax.hlines(50,0,2, color='black', linestyles='dashed', alpha=0.7)

    ax.fill_betweenx([ymin, ymax + 2], 0.75, 1.25, color='lightcoral', alpha=0.5, edgecolor='none')
    ax.set_ylim(ymin, ymax+2)
    ax.set_xlim(0,2)
    
    ax.set_xlabel("z (cm)")
    ax.set_ylabel(r"T "+'('+r'$\circ$'+'C)')

    ax.legend(loc='upper right')
    
    plt.savefig('Euler_implGS_' + str(gamma) + '.png', bbox_inches='tight', dpi=300)

    plt.show()
    
#Comparació i errors

fig_err, ax1 = plt.subplots(figsize=(5,4))
ax1.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')  
ax1.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')


for gamma in gammas:
    filename=f"matriuT_implGS_gamma_{gamma}.txt"
    all_T = np.loadtxt(filename)
    T=all_T[:,-1]


    err = np.abs((T - np.array(T_exacte)) / np.array(T_exacte))

    ymin_err = np.min(err)
    ymax_err = np.max(err)

    ax1.plot(X * L0, err, label=f"$\\gamma = {gamma}$")
    ax1.set_ylim(0, ymax_err+0.000001)

ax1.fill_betweenx([ymin_err, ymax_err+0.000001], 0.75, 1.25, color='lightcoral', alpha=0.5, edgecolor='none')


ax1.set_xlim(0, 2)
ax1.set_xlabel("z (cm)")
ax1.set_ylabel("Error relatiu")
ax1.legend(loc="upper right")
plt.savefig('Erros_implGS_.png', bbox_inches='tight', dpi=300)
plt.show(fig_err)



