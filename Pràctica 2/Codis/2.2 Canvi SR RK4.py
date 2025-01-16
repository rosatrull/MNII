import numpy as np
import matplotlib.pyplot as plt

#IMPORTEM EL VECTOR S(t)
S_t = np.loadtxt("S(t) RK4.txt", delimiter=",", skiprows=1)  #Sense capelra


#PARAMETRES I VARIABLES NECESSÀRIES
G=6.67428e-11                                     
UA=1.496e11                                       
M_S=1.9885e30                                     
T_dia_s=23.9344667*60*60                          
T_any_s=365.25636*24*3600                         
omega=2*np.pi/T_dia_s                             
latitud_c=np.radians(45)                          #latitud de Caldes de Malavella [rad] 41.8356600
longitud_c=np.radians(2.8127200)                  #longitud de Caldes de Malavella [rad]
r_0=UA                                              
t_0=(r_0**3/(M_S*G))**(1/2)
N=int(T_any_s/T_dia_s*24)                         #fem un pas cada hora
dt = T_any_s/(t_0*(N-1))                          #mida del pas temporal
alpha=np.radians(23.44)                           #inclinació eix de la Terra [rad]
M_r_a=np.array([
    [1, 0, 0],
    [0, np.cos(alpha), np.sin(alpha)],
    [0, -np.sin(alpha), np.cos(alpha)]])

#CANVI SISTEMA DE REFERÈNCIA
pos_Caldes_terr=np.zeros((N,3))
a=-1
t=0
for i in S_t:
    t=t+dt
    #Matriu canvi de coordenades
    lamda=latitud_c
    phi_0=longitud_c
    phi=omega*t*t_0+phi_0

    M_r_rad=np.array([
        [-np.cos(phi) * np.sin(lamda), -np.sin(phi) * np.cos(lamda), np.cos(lamda)],
        [-np.sin(phi), np.cos(phi), 0],
        [np.cos(phi) * np.cos(lamda), np.sin(phi) * np.sin(lamda), np.sin(lamda)]
        ])

    #Canvi de base
    a+=1
    pos_Caldes_terr[a]=-np.dot(np.dot(M_r_a,M_r_rad),i)
    print(f"t: {t:.2f}, phi: {phi:.2f}")

R = np.array(pos_Caldes_terr)
filename = f"RK4(t) "
np.savetxt(filename, R, delimiter=",", header="x,y,z", comments="")    

#ORBITA SOL AL VOLTANT DE LA PLACA
x=pos_Caldes_terr[:,0]
y=pos_Caldes_terr[:,1]
z=pos_Caldes_terr[:,2]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z, label="Trajectòria Sol", color='orange')                         # Trajectòries
ax.scatter(0,0,0, color="black", label="placa",s=40)                            # Centre de la Terra
ax.set_title("Posició del Sol resecte la placa")
ax.set_xlabel("Est (m)")
ax.set_ylabel("Nord (m)")
ax.set_zlabel("Radial (m)")
ax.set_zlim(0,1e11)
ax.legend()
plt.show()

#FILTRACIÓ VALORS DE Z
xd=pos_Caldes_terr[:,0]
yd=pos_Caldes_terr[:,1]
zd=pos_Caldes_terr[:,2]
# Filtrar punts on zd > 0
mask = zd > 0  # Això crea una màscara booleana

# Dividir els punts en grups continus
segments = []
current_segment = []

for i in range(len(zd)):
    if mask[i]:  # Si el punt compleix z > 0
        current_segment.append((xd[i], yd[i], zd[i]))
    elif current_segment:  # Si trobem un punt que no compleix voldrà dir que haurem acabat el segment i el guardem a la llista segments
        segments.append(current_segment)
        current_segment = [] #reiniciem la llista per buscar el proxim conjunt de punts (segment) que ho compleix.

# Afegir l'últim segment si n'hi ha
if current_segment:
    segments.append(current_segment)

# Graficar cada segment com una línia independent
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Centre de la Terra (placa)
ax.scatter(0, 0, 0, color="black", label="placa", s=40)

# Traçar cada segment
for segment in segments:
    segment = np.array(segment)  # Convertir a array per facilitar el traçat
    ax.plot(segment[:, 0], segment[:, 1], segment[:, 2], color='orange')

# Configuració del gràfic
ax.set_title("Posició del Sol respecte la placa (z > 0)")
ax.set_xlabel("Est (m)")
ax.set_ylabel("Nord (m)")
ax.set_zlabel("Radial (m)")
ax.set_zlim(0, 1e11)
ax.legend(["placa"])  # Només mostrem la llegenda del punt de la Terra
plt.show()
