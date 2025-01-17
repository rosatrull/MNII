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


#DADES DISTÀNCIA SOL TERRA
d_ST=np.loadtxt(f"RK4.txt", delimiter="," ,skiprows=1)
x_ST=d_ST[:,0]
y_ST=d_ST[:,1]
z_ST=d_ST[:,2]

#CONSTANTS I PARÀMETRES FÍSICS NECESSÀRIS
G=6.67428e-11                                     #cte grav universal [m^3·kg^-1·s^-1]
UA=1.496e11                                       #unitat astronòmica [m] =distància Terra-Sol
M_S=1.9885e30                                     #massa sol [kg]
T_dia_s=23.9344667*60*60                          #dia sideral [s]
T_any_s=365.25636*24*3600                         #any sideral [s]
omega=2*np.pi/T_dia_s                             #freqüència de rotació Terra [rad/s]
R_T=6.371e6                                       #radi de la terra [m]

alpha=np.radians(-23.44)                           #inclinació eix de la Terra [rad]

M_r_a=np.array([
    [np.cos(alpha), 0, np.sin(alpha)],
    [0, 1, 0],
    [-np.sin(alpha), 0, np.cos(alpha)]])          #matriu rotació alpha al voltant d'y (rotem el pla xz)

# línia eix de rotació per després posar-la en el plot
linia = np.linspace(0, 6e6, 100)                   # longitud de la línea
# Coordenadas de la línea
lx = linia * np.sin(alpha)                        # Proyección en el eje x
ly = np.zeros_like(linia)                         # La línea está en el plano xz (y=0)
lz = linia * np.cos(alpha)                        # Proyección en el eje z

#NORMALITZACIÓ
r_0=UA                                              
t_0=(r_0**3/(M_S*G))**(1/2)


#DISCRETITZACIÓ
N=int(T_any_s/T_dia_s*24)                         #fem un pas cada hora
dt = T_any_s/(t_0*(N-1))                          #mida del pas temporal

#CONDICIONS INICIALS
latitud_c=np.radians(41.9831100)                  #latitud de Girona [rad] 41.9831100
longitud_c=np.radians(2.82493)                    #longitud de Girona [rad] 2.82493

#VECTOR S(t)
t=0
p_C_terr=np.zeros((N,3))
pos_Girona=np.zeros((N,3))                        #vector Sol-Girona

# Coordenades esfèriques inicials
theta = np.pi / 2 - latitud_c
phi_0 = longitud_c
R_C_0 = np.array([
    R_T * np.sin(theta) * np.cos(phi_0),
    R_T * np.sin(theta) * np.sin(phi_0),
    R_T * np.cos(theta)
])

# Bucle temporal
for j in range(1, N):
    t = t + dt

    # Longitud actual (amb rotació de la Terra)
    phi = omega * (t*t_0) + phi_0

    # Coordenades cartesianes temporals respecte al centre de la Terra
    x_C = R_T * np.sin(theta) * np.cos(phi)
    y_C = R_T * np.sin(theta) * np.sin(phi)
    z_C = R_T * np.cos(theta)
    R_C_prima = np.array([x_C, y_C, z_C])

    # Rotació per la inclinació de l'eix terrestre
    R_C_0_r = np.dot(M_r_a, R_C_0)  # Posició inicial rotada
    p_C_terr[0] = R_C_0_r

    R_C = np.dot(M_r_a, R_C_prima)  # Evolució temporal rotada
    p_C_terr[j] = R_C
    
    #Vector posició Girona respecte el Sol
    #Posició inicial
    pos_Girona[0][0]=x_ST[0]+p_C_terr[0][0]
    pos_Girona[0][1]=y_ST[0]+p_C_terr[0][1]
    pos_Girona[0][2]=p_C_terr[0][2]
    #Evolució temporal
    pos_Girona[j][0]=x_ST[j]+p_C_terr[j][0]
    pos_Girona[j][1]=y_ST[j]+p_C_terr[j][1]
    pos_Girona[j][2]=z_ST[j]+p_C_terr[j][2]          #recordem que z_ST és zero per tot t ja que haviem col·local el sol a (0,0,0)

#GUARDEM EL VECTOR S(t) EN UN FITXER    
S_t= np.array(pos_Girona)
filename = f"S(t) RK4.txt"
np.savetxt(filename, S_t, delimiter=",", header="x,y,z", comments="")

#ORBITA Girona AL VOLTANT DEL CENTRE DE LA TERRA---------------------------------------------------------------
xr=p_C_terr[:,0]
yr=p_C_terr[:,1]
zr=p_C_terr[:,2]

fig = plt.figure(figsize=(5, 5), dpi = 300)
ax = fig.add_subplot(111, projection='3d')

ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
ax.tick_params(axis='z', which='both', direction='in')

ax.grid(False)

ax.xaxis.pane.set_visible(True)
ax.yaxis.pane.set_visible(True)
ax.zaxis.pane.set_visible(True)  
ax.xaxis.pane.set_edgecolor('grey')
ax.yaxis.pane.set_edgecolor('grey')

ax.xaxis.pane.set_facecolor((0.9, 0.9, 0.9, 1.0))  
ax.yaxis.pane.set_facecolor((0.9, 0.9, 0.9, 1.0))  
ax.zaxis.pane.set_facecolor((0.9, 0.9, 0.9, 1.0))  

ax.scatter(xr[0], yr[0], zr[0], color='red', label="Punt inicial")             # Posició inicial
ax.plot(xr, yr, zr, label="Trajectòria de Girona", color='blue')               # Trajectòries
ax.scatter(0,0,0, color="black", label="Centre Terra",s=40)                    # Centre de la Terra
ax.plot(lx, ly, lz, label='Eix de rotació', color='black', linestyle='--')

#ax.set_title("Posició de Girona respecte el centre de la Terra")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.legend(loc='upper right')
plt.show()


#ORBITA Girona AL VOLTANT DEL SOL--------------------------------------------------------------------------------
xr=pos_Girona[:,0]
yr=pos_Girona[:,1]
zr=pos_Girona[:,2]

fig = plt.figure(figsize=(5, 5), dpi = 300)
ax = fig.add_subplot(111, projection='3d')

ax.tick_params(axis='x', which='both', top=True, labeltop=False, direction='in')
ax.tick_params(axis='y', which='both', right=True, labelright=False, direction='in')
ax.tick_params(axis='z', which='both', direction='in')

ax.grid(False)

ax.xaxis.pane.set_visible(True)
ax.yaxis.pane.set_visible(True)
ax.zaxis.pane.set_visible(True)  
ax.xaxis.pane.set_edgecolor('grey')
ax.yaxis.pane.set_edgecolor('grey')

ax.xaxis.pane.set_facecolor((0.9, 0.9, 0.9, 1.0))  
ax.yaxis.pane.set_facecolor((0.9, 0.9, 0.9, 1.0))  
ax.zaxis.pane.set_facecolor((0.9, 0.9, 0.9, 1.0)) 

ax.scatter(xr[0], yr[0], zr[0], color='red', label="Punt inicial (x, y, z)")   # Posició inicial
ax.plot(xr, yr, zr, label="Trajectòria (x, y, z)", color='blue')               # Trajectòries
ax.scatter(0,0,0, color="orange", label="Sol",s=40)                            # Sol

#ax.set_title("Posició de Girona respecte el Sol")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_zlim(-1.5e11,1.5e11)
ax.legend()
plt.show()
