import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
size = 100

# température avec échange par convection avec l'extérieur
Tinf = 300

# delta x correspondant à la longueur de controle
dx = .5

# conductivité thermique du matériau
k = 80

# puissance volumique
q = 100

# coefficient de transfert par convection
h = 10

def solve():
    t1=time.perf_counter()
    ax=Axes3D(plt.figure())
    A = np.zeros((size**2, size**2))
    B = np.zeros((size**2))
    
    for k in range(0, size**2):
        
        # coordonnées du point
        x = k % size
        y = k // size
        
        # nature du point
        # point = findNature(x, y)
        
        # A[k] = point.getEquation()
        # B[k] = -q * dx**2 + point.getConvectionTerm()
    
    C=np.linalg.solve(A, B)
    print(C)
    E=np.zeros((size,size))
    for i in range(size) :
        for j in range(size) :
            E[i][j]=C[i+size*j]
    X,Y=np.mgrid[0:size:1,0:size:1]
    ax.plot_wireframe(X,Y,E)
    t2=time.perf_counter()
    print(t2-t1)
    plt.show()
    
# donne le numero de coefficient dans l'équation linéaire, en fonction de x, y, et z
# x, y, z concaténés représentent l'écriture de getIndex(x, y, z) en base 'size'
def getIndex(x, y, z):
    return x + y * size + z * size**2

class Point:
    def __init__(self, x, y, z, c):
        
        # coordonées du point
        self.x = x
        self.y = y
        self.z = z
        
        # contacts avec les blocs élémentaires adjacents
        # liste d'offsets, len(c) = nb de contacts
        # (6 - len(c)) = nb de contacts convectifs
        self.c = c
    
    def getEquation(self):
        A = np.zeros(size**3)
        index = getIndex(self.x, self.y, self.z)
        
        # la température du bloc intervient autant de fois dans l'équation linéaire
        # qu'il y a de contacts avec les blocs adjacents, d'où len(c) * k
        
        # (6 - len(c)) * h * dc représente le fait que la temp du bloc intervient
        # dans le delta de température dans le calcul du flux thermique convectif,
        # d'où (6 - len(c)) * h * dx
        A[index] = -len(self.c) * k - (6 - len(self.c)) * h * dx
        
        # appliquer les coefficients dans l'équation linéaire car les températures des
        # blocs adjacents interviennent avec un facteur k dans le transfert par conduction
        for a in self.c:
            A[index + a] = k
    
    # terme de convection à droite de l'égalité, il ne faut pas oublier
    # de lui sommer le terme de puissance volumique.
    def getConvectionTerm(self):
        return - (6 - len(self.c)) * h * dx * Tinf

# Un centre est défini par aucun contact convectif
# tous les contacts sont conductifs et il n'y a aucune ambiguité sur l'orientation des contacts
class Centre(Point):
    
    # Aucun contact convectif
    def __init__(self, x, y, z):
        super().__init__(x, y, z, [-1, 1, size, -size, size**2, -size**2])

# Un "coin" possède 3 contacts convectifs
# Un "bord" possède 2 contacts convectifs
# Une "surface" possède 3 contacts convectifs
# Une "lame" possède 4 contacts convectifs


# Un programme tiers codé sur Minecraft permet d'analyser une construction, et de former
# une array de chiffres qui correspondent à des types de points, qui sont ensuite lus
# par le logiciel Python pour former la géométrie