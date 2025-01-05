import numpy as np
import matplotlib.pyplot as plt

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

#CONSTANTS I PARÀMETRES FÍSICS NECESSÀRIS-----------------------------------------------------------------------------
G=6.67428e-11                                     #cte grav universal [m^3·kg^-1·s^-1]
UA=1.496e11                                       #unitat astronòmica [m] =distància Terra-Sol
M_S=1.9885e30                                     #massa sol [kg]
T_dia_s=23.9344667*60*60                          #dia sideral [s]
T_any_s=365.25636*24*3600                         #any sideral [s]
omega=2*np.pi/T_dia_s                             #freqüència de rotació Terra [rad/s]

#NORMALITZACIÓ
r_0=UA                                              
t_0=(r_0**3/(M_S*G))**(1/2)


#DISCRETITZACIÓ
N=int(T_any_s/T_dia_s*24)                         #fem un pas cada hora
dt = T_any_s/(t_0*(N-1))                          #mida del pas temporal

#CONDICIONS INICIALS------------------------------------------------------------------------------------------------------
d_af=1.017*UA/r_0                                 #distància Terra-Sol a l'afeli [m]
v_af=-(29300*t_0)/r_0                             #velocitat Terra a l'afeli [m/s] !Component y!

pos_f=np.zeros((N,3))                             #(x,y,z) a t_n+1
pos_f[0,0]=d_af
pos_f[0,1]=0

vel_f=np.zeros((N,3))                             #(vx,vy,vz) a t_n+1
vel_f[0,0]=0
vel_f[0,1]=v_af

t=0
#MÈTODE D'EULER----------------------------------------------------------------------------------------------------------
def acceleracio_x(x,y):
    ax = - x / (x**2 + y**2)**(3/2)
    return ax
def acceleracio_y(x,y):
    ay = - y / (x**2 + y**2)**(3/2)
    return ay

for j in range(1,N):

    #Actualitzem velocitats
    vel_f[j,0] = vel_f[j-1,0] + acceleracio_x(pos_f[j-1,0], pos_f[j-1,1]) * dt
    vel_f[j,1] = vel_f[j-1,1] + acceleracio_y(pos_f[j-1,0], pos_f[j-1,1]) * dt
    #Actualitzem posicions
    pos_f[j,0] = pos_f[j-1,0] + vel_f[j-1,0] * dt 
    pos_f[j,1] = pos_f[j-1,1] + vel_f[j-1,1] * dt 
   

#GRAFIQUEM ORBITA TERRA AL VOLTANT DEL SOL
xd=pos_f[:, 0]*r_0
yd=pos_f[:, 1]*r_0
#zd=pos_f[:, 2]*r_0

fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
ax.plot(xd, yd, label='Orbita de la Terra al voltant del Sol amb Euler', color='blue')

ax.scatter(0, 0, s=40, label="Sol", color="orange")

#ax.set_ylim(min(yd)-1e11, max(yd)+1e11)
#ax.set_xlim(min(xd)-1e11, max(xd)+1e11)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.legend(loc='upper right')

# Guardar el gràfic
plt.savefig('Orbita_Terra_Euler', bbox_inches='tight', dpi=300)
plt.show()
