import numpy as np
L0=2

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


gammas = [0.25] 

n = 101  #Nombre de punts d'espai
X = np.linspace(0, 1, n)  # Mallat espacial (adimensional)
dx = 1 / (n - 1)  # Amplada del mallat espacial
t_a = 0.025  #

for gamma in gammas:
    filename=f"matriuT_expl_gamma_{gamma}.txt"
    T = np.loadtxt(filename)
    dt = gamma * dx**2
    m = int(t_a / dt) + 1  
    ymin=np.min(T)-273.15
    ymax=np.max(T)-273.15
    print(m)

    num_images = int(m/6)
    
    indices = np.linspace(0, m-1, num_images, dtype=int)


    #volem fer una evolució temporal, pel que a cada plot cal plotejar cada columna. 
    fig, ax = plt.subplots(figsize=(5,4))
    ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')  
    ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')  

    ax.set_ylim(ymin,ymax+1)
    ax.set_xlabel("z (cm)")
    ax.set_ylabel(r"T "+'('+r'$\circ$'+'C)')
    ax.set_title(f"Mètode d'Euler Explícit. Gamma={gamma}")
    ax.legend_ = None 

    ax.fill_betweenx([ymin, ymax + 1], 0.75, 1.25, color='lightcoral', alpha=0.5, edgecolor='none')
    # Loop over each time step and add the new curve while updating the legend
    for k in indices:
        # Plot the current time step curve
        current_line, = ax.plot(
            X * L0, T[:, k]-273.15,
            label=f"$\\hat{{t}}={k*dt:.4f}$",
            alpha=((0.9/m)*k+0.1),
            color="blue"
        )
        
        # Remove the old legend and add a new one with only the current curve
        ax.legend([current_line], [f"$\\hat{{t}}={k*dt:.4f}$"], loc='upper right', fontsize=10)

        # Save the plot for this time step
        plt.savefig(f"temperature_plot_gamma_{gamma}_t_{k*dt:.4f}.png", dpi=150)

    # Close the figure after saving all images
    plt.close(fig)

