import numpy as np
import matplotlib.pyplot as plt

#CONSTANTS I PARÀMETRES FÍSICS NECESSÀRIS--------------------------------------------------------
G=6.67428e-11                                     #cte grav universal [m^3·kg^-1·s^-1]
UA=1.496e11                                       #unitat astronòmica [m] =distància Terra-Sol
M_S=1.9885e30                                     #massa sol [kg]
#M_T=5.97*(10**(24))                              #massa terra [kg]
R_T=6.371e6                                       #radi de la terra [m]

alpha=np.radians(23.44)                           #inclinació eix de la Terra [rad]
M_r_a=np.array([
    [1, 0, 0],
    [0, np.cos(alpha), np.sin(alpha)],
    [0, -np.sin(alpha), np.cos(alpha)]])          #matriu rotació alpha al voltant d'x

T_dia_s=23.9344667*60*60                          #dia sideral [s]
T_any_s=365.25636*T_dia_s                         #any sideral [s]
omega=2*np.pi/T_dia_s                             #freqüència de rotació Terra [rad/s]

#DISCRETITZACIÓ----------------------------------------------------------------------------------
N=365*24                                          #fem un pas cada hora
dt=1/(N-1)                                        #mida del pas temporal

#NORMALITZACIÓ-----------------------------------------------------------------------------------
t_0=365*24*60*60                                  #Normalització temps
r_0=(G*M_S*t_0**2)**(1/3)

#CONDICIONS INICIALS-----------------------------------------------------------------------------
d_af=1.017*UA/r_0                                 #distància Terra-Sol a l'afeli [m]
v_af=-(29300*t_0)/r_0                             #velocitat Terra a l'afeli [m/s] !Component y!
latitud_c=np.radians(41.8356600)                  #latitud de Caldes de Malavella [rad] 
longitud_c=np.radians(2.8127200)                  #longitud de Caldes de Malavella [rad]

pos_i=np.zeros((N,3))                             #(x,y,z) a t n
pos_f=np.zeros((N,3))                             #(x,y,z) a t n+1
pos_f[0,0]=d_af
pos_f[0,1]=0

vel_i=np.zeros((N,3))                             #(vx,vy,vz) a t n
vel_f=np.zeros((N,3))                             #(vx,vy,vz) a t n
vel_f[0,0]=0
vel_f[0,1]=v_af

t=0

p_C_terr=np.zeros((N,3))
pos_Caldes=np.zeros((N,3))


#RUNGE KUTTA-------------------------------------------------------------------------------------
# RUNGE-KUTTA
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
    
    #Vector posició Caldes respect el centre de la Terra
    theta=np.pi/2-latitud_c
    phi_0=longitud_c
    phi=omega*t*t_0+phi_0
    
    #Posició inical
    R_C_0=[R_T*np.cos(theta)*np.cos(phi_0),
                 R_T*np.cos(theta)*np.sin(phi_0),
                 R_T*np.sin(theta)]
    
    #Evolució temporal
    x_C=R_T*np.cos(theta)*np.cos(phi)
    y_C=R_T*np.cos(theta)*np.sin(phi)
    z_C=R_T*np.sin(theta)
    R_C_prima=np.array((x_C,y_C,z_C))
    
    #Rotació eix Terra
    R_C_0_r=np.dot(M_r_a,R_C_0)
    p_C_terr[0]=R_C_0_r
    R_C=np.dot(M_r_a,R_C_prima)
    p_C_terr[j]=R_C
    
    #Vector posició Caldes respect Sol
    #Posició inicial
    pos_Caldes[0][0]=pos_f[0][0]*r_0+R_T*np.cos(theta)*np.cos(phi_0)
    pos_Caldes[0][1]=pos_f[0][1]*r_0+R_T*np.cos(theta)*np.sin(phi_0)
    pos_Caldes[0][2]=z_C
    #Posicions següents
    pos_Caldes[j][0]=pos_f[j][0]*r_0+x_C
    pos_Caldes[j][1]=pos_f[j][1]*r_0+y_C
    pos_Caldes[j][2]=z_C


#ORBITA TERRA AL VOLTANT DEL SOL---------------------------------------------------------------    
xd=pos_f[:, 0]*r_0
yd=pos_f[:, 1]*r_0
zd=pos_f[:, 2]

plt.figure(figsize=(8, 8))
plt.plot(xd,yd,label="orbita Terra")
plt.scatter(0,0, color="orange",s=40)
plt.legend()
plt.show()

#ORBITA CALDES AL VOLTANT DEL CENTRE DE LA TERRA--------------------------------------------------------------
xr=p_C_terr[:,0]
yr=p_C_terr[:,1]
zr=p_C_terr[:,2]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xr[0], yr[0], zr[0], color='red', label="Punt inicial (x, y, z)")   # Posició inicial
ax.plot(xr, yr, zr, label="Trajectòria (x, y, z)", color='blue')               # Trajectòries
ax.scatter(0,0,0, color="black", label="Centre Terra",s=40)                            # Centre de la Terra
ax.set_title("Posició de Caldes respecte el centre de la Terra")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.legend()
plt.show()

#ORBITA CALDES AL VOLTANT DEL SOL--------------------------------------------------------------
xr=pos_Caldes[:,0]
yr=pos_Caldes[:,1]
zr=pos_Caldes[:,2]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xr[0], yr[0], zr[0], color='red', label="Punt inicial (x, y, z)")   # Posició inicial
ax.plot(xr, yr, zr, label="Trajectòria (x, y, z)", color='blue')               # Trajectòries
ax.scatter(0,0,0, color="orange", label="Sol",s=40)                            # Centre de la Terra
ax.set_title("Posició de Caldes respecte el Sol")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_zlim(-1.5e11,1.5e11)
ax.legend()
plt.show()



#CANVI SISTEMA DE REFERÈNCIA------------------------------------------------------------------
pos_Caldes_terr=np.zeros((N,3))
a=-1
t=0
for i in pos_Caldes:
    t=t+dt
    #Matriu canvi de coordenades
    theta=np.pi/2-latitud_c 
    phi=omega*t*t_0+phi_0
    #M_r_rad=np.array([
    #    [-np.sin(phi), -np.sin(theta) * np.cos(phi), np.cos(theta) * np.cos(phi)],
    #    [np.cos(phi), -np.sin(theta) * np.sin(phi), np.cos(theta) * np.sin(phi)],
    #    [0,np.cos(theta), np.sin(theta)]
    #    ])
    M_r_rad = np.array([
    [np.cos(phi), -np.sin(phi), 0],
    [np.sin(phi) * np.cos(theta), np.cos(phi) * np.cos(theta), -np.sin(theta)],
    [np.sin(phi) * np.sin(theta), np.cos(phi) * np.sin(theta), np.cos(theta)]
    ])
    #Canvi de base
    a+=1
    pos_Caldes_terr[a]=-np.dot(M_r_rad,i)
    
#ORBITA SOL AL VOLTANT DE LA PLACA--------------------------------------------------------------
x=pos_Caldes_terr[:,0]
y=pos_Caldes_terr[:,1]
z=pos_Caldes_terr[:,2]

dia=pos_Caldes_terr[:25]
xd=dia[:,0]
yd=dia[:,1]
zd=dia[:,2]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(y, x, z, label="Trajectòria Sol", color='orange')                       # Trajectòries
ax.scatter(0,0,0, color="black", label="placa",s=40)                            # Centre de la Terra
plt.plot(yd,xd,zd, color='blue', label='1r dia')
ax.set_title("Posició del Sol resecte la placa")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.legend()
plt.show()

#PROJECCIONS 2D
plt.plot(x,y,label="x,y")
plt.scatter(0,0, color="orange",s=40)
plt.legend()
plt.show()

plt.plot(y,z, label="y,z")
plt.scatter(0,0, color="orange",s=40)
plt.legend()
plt.show()

plt.plot(x,z,label="x,y")
plt.scatter(0,0, color="orange",s=40)
plt.legend()
plt.show()

#DIES POSSIBLES----------------------------------------------------------------------------------
b=0 
S=[]
S_tot=[]
for i in pos_Caldes_terr:
    b+=1
    S.append(i)
    if b==25:
        S_tot.append(S)
        S=[]
        b=0

#4 juliol (distància afeli)
q_juliol=np.array(S_tot[0])
x1= q_juliol[:,0]
y1= q_juliol[:,1]
z1= q_juliol[:,2]

stjoan=np.array(S_tot[150])
x2= stjoan[:,0]
y2= stjoan[:,1]
z2= stjoan[:,2]

nitbona=np.array(S_tot[336])
x3= nitbona[:,0]
y3= nitbona[:,1]
z3= nitbona[:,2]

dia_r=np.array(S_tot[100])
x4= dia_r[:,0]
y4= dia_r[:,1]
z4= dia_r[:,2]

#Graficar 3D
fig = plt.figure(figsize=(13, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot(y1,x1,z1, color="orange",label="1 gener")
ax.plot(y2,x2,z2, color="c",label="24 juny")
#ax.plot(y3,x3,z3, color=(128/255, 0/255, 0/255), label="24 desembre")
#ax.plot(y4,x4,z4, color="green", label="dia random")
ax.scatter(0,0,0, s=30, color="black", label="Placa")
# Etiquetes i títol
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("radial (m)")
ax.set_title("Sol des de la Placa")
ax.legend()
plt.show()
