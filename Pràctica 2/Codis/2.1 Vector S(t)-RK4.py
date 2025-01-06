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
d_ST=np.loadtxt(f"RK4.txt",delimiter=",",skiprows=1)
x_ST=d_ST[:,0]
y_ST=d_ST[:,1]

#CONSTANTS I PARÀMETRES FÍSICS NECESSÀRIS
G=6.67428e-11                                     #cte grav universal [m^3·kg^-1·s^-1]
UA=1.496e11                                       #unitat astronòmica [m] =distància Terra-Sol
M_S=1.9885e30                                     #massa sol [kg]
T_dia_s=23.9344667*60*60                          #dia sideral [s]
T_any_s=365.25636*24*3600                         #any sideral [s]
omega=2*np.pi/T_dia_s                             #freqüència de rotació Terra [rad/s]
R_T=6.371e6                                       #radi de la terra [m]

alpha=np.radians(23.44)                           #inclinació eix de la Terra [rad]
M_r_a=np.array([
    [1, 0, 0],
    [0, np.cos(alpha), np.sin(alpha)],
    [0, -np.sin(alpha), np.cos(alpha)]])          #matriu rotació alpha al voltant d'x

# línia eix de rotació per després posar-la en el plot
linia = np.linspace(0, 6e6, 100)                   # longitud de la línea
# Coordenadas de la línea
lx = np.zeros_like(linia)                         # La línea está en el plano yz (x=0)
ly = linia * np.sin(alpha)                        # Proyección en el eje y
lz = linia * np.cos(alpha)                        # Proyección en el eje z

#NORMALITZACIÓ
r_0=UA                                              
t_0=(r_0**3/(M_S*G))**(1/2)


#DISCRETITZACIÓ
N=int(T_any_s/T_dia_s*24)                         #fem un pas cada hora
dt = T_any_s/(t_0*(N-1))                          #mida del pas temporal

#CONDICIONS INICIALS
latitud_c=np.radians(45)                          #latitud de Caldes de Malavella [rad] 41.8356600
longitud_c=np.radians(2.8127200)                  #longitud de Caldes de Malavella [rad]


#VECTOR S(t)
t=0
p_C_terr=np.zeros((N,3))
pos_Caldes=np.zeros((N,3))
for j in range(1, N):
    t = t + dt

    #Vector posició Caldes respect el centre de la Terra
    theta=np.pi/2-latitud_c
    phi_0=longitud_c
    phi=omega*t*t_0+phi_0
    
    #Posició inical
    R_C_0=[R_T*np.sin(theta)*np.cos(phi_0), R_T*np.sin(theta)*np.sin(phi_0), R_T*np.cos(theta)]
    #Evolució temporal
    x_C=R_T*np.sin(theta)*np.cos(phi)
    y_C=R_T*np.sin(theta)*np.sin(phi)
    z_C=R_T*np.cos(theta)
    R_C_prima=np.array((x_C,y_C,z_C))
    
    #Rotació eix Terra
    #Posició inicial
    R_C_0_r=np.dot(M_r_a,R_C_0)
    p_C_terr[0]=R_C_0_r
    #Evolució temporal
    R_C=np.dot(M_r_a,R_C_prima)
    p_C_terr[j]=R_C
    
    #Vector posició Caldes respect Sol
    #Posició inicial
    pos_Caldes[0][0]=x_ST[0]+p_C_terr[0][0]
    pos_Caldes[0][1]=y_ST[0]+p_C_terr[0][1]
    pos_Caldes[0][2]=p_C_terr[0][2]
    #Evolució temporal
    pos_Caldes[j][0]=x_ST[j]+p_C_terr[j][0]
    pos_Caldes[j][1]=y_ST[j]+p_C_terr[j][1]
    pos_Caldes[j][2]=p_C_terr[j][2]
 
#GUARDEM EL VECTOR S(t) EN UN FITXER    
S_t= np.array(pos_Caldes)
filename = f"S(t) RK4.txt"
np.savetxt(filename, S_t, delimiter=",", header="x,y,z", comments="")
      
#ORBITA CALDES AL VOLTANT DEL CENTRE DE LA TERRA---------------------------------------------------------------
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
ax.plot(xr, yr, zr, label="Trajectòria de Caldes", color='blue')               # Trajectòries
ax.scatter(0,0,0, color="black", label="Centre Terra",s=40)                    # Centre de la Terra
ax.plot(lx, ly, lz, label='Eix de rotació', color='black', linestyle='--')

ax.set_title("Posició de Caldes respecte el centre de la Terra")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.legend(loc='upper right')
plt.show()


#ORBITA CALDES AL VOLTANT DEL SOL--------------------------------------------------------------------------------
xr=pos_Caldes[:,0]
yr=pos_Caldes[:,1]
zr=pos_Caldes[:,2]

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

ax.set_title("Posició de Caldes respecte el Sol")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_zlim(-1.5e11,1.5e11)
ax.legend()
plt.show()
