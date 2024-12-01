import numpy as np
import matplotlib.pyplot as plt

# Paràmetres inicials
n = 101  
X = np.linspace(0, 1, n) 
dx = 1 / (n - 1)  
t_a = 0.025  
gammas = [0.5, 1]  
T_0 = (0.02**2 * 944000) / (0.56)
solucions = []

# Iteració per gamma
for gamma in gammas:
    dt = gamma * dx**2
    m = int(t_a / dt) + 1                              
    T = np.zeros((n, m))                                
    T[:, 0] = 309.65 / T_0                              
    T[0, :] = 309.65 / T_0                              
    T[-1, :] = 309.65 / T_0                             
    
    for i in range(1, T.shape[1]):                      
        T_prev = T[:, i-1]                              
        T_nou_i = T_prev.copy()                         
        T_nou = T_prev.copy()                           

        error = 1e6                                     
        while error > 1e-6:                             
            for j in range(T.shape[0]):                 
                if j == 0 or j == n - 1:
                    T_nou[j] = 309.65 / T_0
                else:
                    T_nou[j] = (T_prev[j] + gamma * (T_nou_i[j + 1] + T_nou_i[j - 1]) + dt) / (1 + 2 * gamma)
    
            error = np.linalg.norm(T_nou - T_nou_i, ord=np.inf)  
            if error > 1e-6:
                T_nou_i = T_nou.copy()  
    
        T[:, i] = T_nou
    solucions.append((T[:, -1] * T_0, gamma))              


def T_tilda(t, z, T_c, N_terms=100):
    coef = 4 / (np.pi**3)
    summation = 0
    for n in range(N_terms):
        k = 2 * n + 1  
        term = (1 - np.exp(-(k**2) * (np.pi**2) * t)) / (k**3) * np.sin(k * np.pi * z)
        summation += term
    return T_c + coef * summation

# Calculem la solució analítica
t = 0.025  
z = np.linspace(0, 1, n)  
T_c = 309.65 / T_0  
N_terms = 1000  
T_exacte = [T_tilda(t, zi, T_c) * T_0 for zi in z]

# Error absolut per gamma
fig, ax1 = plt.subplots(figsize=(5, 4), dpi=500)
ax1.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')  
ax1.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')

for i, gamma in enumerate(gammas):
    temperatures = solucions[i][0]
    err_abs = np.abs((temperatures - T_exacte))
    ax1.plot(X*2, err_abs, label=f"$\\gamma = {gamma}$")


ax1.set_xlabel("z (cm)")
ax1.set_ylabel("Error absolut")
ax1.legend(loc="upper right")
ax1.set_xlim(0, 2)
ax1.set_ylim(min(err_abs),max(err_abs)+0.01)
ax1.fill_betweenx([min(err_abs),max(err_abs)+0.01], 0.75 , 1.25 , 
                     color='lightcoral', alpha=0.5, edgecolor='none')
plt.savefig('Erros_implCncpt_.png', bbox_inches='tight', dpi=300)


plt.show()
for i, gamma in enumerate(gammas):
    temperatures = np.array(solucions[i][0]) - 273.15  
    T_exacte_C = np.array(T_exacte) - 273.15  
    fig, ax = plt.subplots(figsize=(5, 4))
    
    ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')  
    ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
    
    ax.fill_betweenx([np.min(temperatures), np.max(temperatures)+2], 0.75 , 1.25 , 
                     color='lightcoral', alpha=0.5, edgecolor='none')
    
    ax.plot(X * 2, T_exacte_C, label='Solució analítica')
    ax.scatter(X * 2, temperatures, label=f"$\\gamma = {gamma}$\n$\\hat{{t}}=0.025$")
    ax.hlines(50,0,2, color='black', linestyles='dashed', alpha=0.7)
    ax.set_xlim(0,2)
    ax.set_ylim(np.min(temperatures), np.max(temperatures)+2)
    ax.set_xlabel("z (cm)")
    ax.set_ylabel(r"T "+'('+r'$\circ$'+'C)')
    ax.legend(loc="upper right")
    plt.savefig('Euler_implConcpt_' + str(gamma) + '.png', bbox_inches='tight', dpi=300)
    
    plt.show()
