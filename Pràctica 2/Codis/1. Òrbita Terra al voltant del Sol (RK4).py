import numpy as np
import matplotlib.pyplot as plt

#CONSTANTS I PARÀMETRES FÍSICS NECESSÀRIS
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

#CONDICIONS INICIALS
d_af=1.017*UA/r_0                                 #distància Terra-Sol a l'afeli [m]
v_af=-(29300*t_0)/r_0                             #velocitat Terra a l'afeli [m/s] !Component y!

pos_i=np.zeros((N,3))                             #(x,y,z) a t_n
pos_f=np.zeros((N,3))                             #(x,y,z) a t_n+1
pos_f[0,0]=d_af
pos_f[0,1]=0

vel_i=np.zeros((N,3))                             #(vx,vy,vz) a t_n
vel_f=np.zeros((N,3))                             #(vx,vy,vz) a t_n+1
vel_f[0,0]=0
vel_f[0,1]=v_af


#RUNGE-KUTTA 4
t=0
p_C_terr=np.zeros((N,3))
pos_Caldes=np.zeros((N,3))

for j in range(1, N):
    t = t + dt

    # Actualitzem k1
    k1x = dt * vel_f[j-1][0]
    k1y = dt * vel_f[j-1][1]
    k1vx = -dt * pos_f[j-1][0] / (pos_f[j-1][0]**2 + pos_f[j-1][1]**2)**(3/2)
    k1vy = -dt * pos_f[j-1][1] / (pos_f[j-1][0]**2 + pos_f[j-1][1]**2)**(3/2)

    # Actualitzem k2
    k2x = dt * (vel_f[j-1][0] + k1vx / 2)
    k2y = dt * (vel_f[j-1][1] + k1vy / 2)
    k2vx = -dt * (pos_f[j-1][0] + k1x / 2) / (((pos_f[j-1][0] + k1x / 2)**2 + (pos_f[j-1][1] + k1y / 2)**2)**(3/2))
    k2vy = -dt * (pos_f[j-1][1] + k1y / 2) / (((pos_f[j-1][0] + k1x / 2)**2 + (pos_f[j-1][1] + k1y / 2)**2)**(3/2))

    # Actualitzem k3
    k3x = dt * (vel_f[j-1][0] + k2vx / 2)
    k3y = dt * (vel_f[j-1][1] + k2vy / 2)
    k3vx = -dt * (pos_f[j-1][0] + k2x / 2) / (((pos_f[j-1][0] + k2x / 2)**2 + (pos_f[j-1][1] + k2y / 2)**2)**(3/2))
    k3vy = -dt * (pos_f[j-1][1] + k2y / 2) / (((pos_f[j-1][0] + k2x / 2)**2 + (pos_f[j-1][1] + k2y / 2)**2)**(3/2))

    # Actualitzem k4
    k4x = dt * (vel_f[j-1][0] + k3vx)
    k4y = dt * (vel_f[j-1][1] + k3vy)
    k4vx = -dt * (pos_f[j-1][0] + k3x) / (((pos_f[j-1][0] + k3x)**2 + (pos_f[j-1][1] + k3y)**2)**(3/2))
    k4vy = -dt * (pos_f[j-1][1] + k3y) / (((pos_f[j-1][0] + k3x)**2 + (pos_f[j-1][1] + k3y)**2)**(3/2))

    # Posició terra respecte el Sol
    pos_f[j][0] = pos_f[j-1][0] + (k1x + 2*k2x + 2*k3x + k4x) / 6
    pos_f[j][1] = pos_f[j-1][1] + (k1y + 2*k2y + 2*k3y + k4y) / 6
    vel_f[j][0] = vel_f[j-1][0] + (k1vx + 2*k2vx + 2*k3vx + k4vx) / 6
    vel_f[j][1] = vel_f[j-1][1] + (k1vy + 2*k2vy + 2*k3vy + k4vy) / 6 
    
    
#GRAFIQUEM ORBITA TERRA AL VOLTANT DEL SOL
xd=pos_f[:, 0]*r_0
yd=pos_f[:, 1]*r_0
zd=pos_f[:, 2]*r_0

plt.figure(figsize=(8, 8))
plt.plot(xd,yd,label="orbita Terra")
plt.scatter(0,0, color="orange",s=40)
plt.legend()
plt.show()