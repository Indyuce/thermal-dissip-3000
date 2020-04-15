import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# conductivité thermique de l'aluminium à 20°C (W/mK)
conduc = 237

# paramètres géométriques
e = 0.01 # epaisseur d'un bloc (m)
S = e**2 # surface en contact entre deux blocs (m²)

# masse totale du dissipateur (kg)
m = 0.3

# résistance thermique d'un bloc
R = e / (conduc * S)

# capacité thermique massique (J/kg/K)
# permet d'établir le lien entre energie thermique et température
c = 897

# paramètres de la simulation
time = 300
iterations = 1000
dt = time / iterations
size = 20

# température initiale de la source au milieu
tempInit = 80

# constante de simulation pour limiter les calculs
cste = dt / (R * m * c)

def simulation():
    
    # espace à deux dimensions
    T = np.zeros((size, size))

    # source de chaleur initiale
    T[size // 2][size // 2] = tempInit

    for k in range(iterations):

        # la source de chaleur emet une certaine energie thermique chaque seconde
        # T[size // 2][size // 2] += temp

        B = np.array(T)

        for k in range(2, size - 1):
            for j in range(2, size - 1):
                temp = T[k][j]
                B[k][j] += (T[k + 1][j] + T[k - 1][j] + T[k][j + 1] + T[k][j - 1] - 4 * temp) * cste

        T = B

    X, Y = np.mgrid[0:size:1, 0:size:1]
    Axes3D(plt.figure()).plot_wireframe(X, Y, T)    
    plt.show()

simulation()