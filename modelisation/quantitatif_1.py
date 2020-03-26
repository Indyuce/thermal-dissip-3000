import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# conductivité thermique de l'aluminium à 20°C (W/mK)
conduc = 185

# paramètres géométriques
e = 0.01 # epaisseur d'un bloc (m)
S = e**2 # surface en contact entre deux blocs (m²)

# masse totale du dissipateur (kg)
m = 0.3

# résistance thermique d'un bloc
R = e / (conduc * S)

# capacité thermique (J/kg/K)
# permet d'établir le lien entre energie thermique et température
c = 900

# paramètres de la simulation
time = 100
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
                B[k][j] += (T[k + 1][j] - T[k][j]) * cste
                B[k][j] += (T[k - 1][j] - T[k][j]) * cste
                B[k][j] += (T[k][j + 1] - T[k][j]) * cste
                B[k][j] += (T[k][j - 1] - T[k][j]) * cste

        T = B

    # conservation de l'énergie thermique
    s = 0
    for k in range(0, size):
        for l in range(0, size):
            s += T[k][l]
    
    if (tempInit - s) / tempInit > 0.01:
        print('Conservation de l''énergie non vérifiée:', tempInit, s)

    X, Y = np.mgrid[0:size:1, 0:size:1]
    Axes3D(plt.figure()).plot_wireframe(X, Y, T)
    plt.show()

simulation()