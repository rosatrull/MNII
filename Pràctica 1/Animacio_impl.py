####Animacio

import numpy as np
L0=2

#intentem no tocar això, així tots els plots quedaràn amb mides iguals. 

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer modern"],  
    "axes.labelsize": 14,     
    "axes.linewidth": 0.5,    
    "xtick.labelsize": 12,    
    "ytick.labelsize": 12,    
    "legend.fontsize": 10,    
    "legend.handlelength": 1.0,
    "lines.linewidth": 1,    
    "lines.markersize": 3     
})


gammas = [0.5] 

n = 101  #Nombre de punts d'espai
X = np.linspace(0, 1, n)  # Mallat espacial (adimensional)
dx = 1 / (n - 1)  # Amplada del mallat espacial
t_a = 0.025  #

for gamma in gammas:
    filename=f"matriuT_implGS_gamma_{gamma}.txt"
    T = np.loadtxt(filename)
    dt = gamma * dx**2
    m = int(t_a / dt) + 1  
    ymin=np.min(T)
    ymax=np.max(T)
    print(m)

    num_images = 100
    
    indices = np.linspace(0, m-1, num_images, dtype=int)


    #volem fer una evolució temporal, pel que a cada plot cal plotejar cada columna. 
    fig, ax = plt.subplots(figsize=(5,4))
    ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')  
    ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')  

    ax.set_ylim(ymin,ymax+1)
    ax.set_xlabel("z (cm)")
    ax.set_ylabel(r"T "+'('+r'$\circ$'+'C)')
    ax.set_title(f"Mètode d'Euler Implicit. Gamma={gamma}")
    ax.legend_ = None 

    ax.fill_betweenx([ymin, ymax + 1], 0.75, 1.25, color='lightcoral', alpha=0.5, edgecolor='none')

    for k in indices:

        current_line, = ax.plot(
            X * L0, T[:, k],
            label=f"$\\hat{{t}}={k*dt:.4f}$",
            alpha=((0.9/m)*k+0.1),
            color="indigo"
        )
        

        ax.legend([current_line], [f"$\\hat{{t}}={k*dt:.4f}$"], loc='upper right', fontsize=10)

  
        plt.savefig(f"temperature_plot_gamma_{gamma}_t_{k*dt:.4f}.png", dpi=150)

    
    plt.close(fig)


