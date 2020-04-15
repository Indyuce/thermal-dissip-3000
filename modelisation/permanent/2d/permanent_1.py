import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from abc import ABC, abstractmethod

size = 10

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

class Point(ABC):
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    @abstractmethod
    def getEquation(self):
        pass
    
    # terme de convection à droite de l'égalité
    @abstractmethod
    def getConvectionTerm(self):
        pass

class Centre(Point):
    def __init__(self, x, y):
        super().__init__(x, y)
        
    def getEquation(self):
        A = np.zeros(size**2)
        index = size * self.x + self.y
        A[index] = -4 * k
        A[index - 1] = k
        A[index + 1] = k
        A[index - size] = k
        A[index + size] = k
        return A
    
    def getTerme(self):
        return 0

## coins

# coté et hauteur de métal en contact
# côté = 1 pour gauche, -1 pour droite
# hauteur = 1 pour supérieur, -1 pour inférieur
class Coin(Point):
    def __init__(self, x, y, cote, hauteur):
        super().__init__(x, y)
        self.cote = cote
        self.hauteur = hauteur
        
    def getEquation(self):
        A = np.zeros(size**2)
        index = size * self.x + self.y
        A[index + self.cote] = k
        A[index - size * self.hauteur] = k
        A[index] = -2 * k - 2 * h * dx
        return A
    
    def getConvectionTerm(self):
        return -2 * h * dx * Tinf

class CoinInferieurGauche(Coin):
    def __init__(self, x, y):
        super().__init__(x, y, 1, 1)

class CoinInferieurDroit(Coin):
    def __init__(self, x, y):
        super().__init__(x, y, -1, 1)

class CoinSuperieurGauche(Coin):
    def __init__(self, x, y):
        super().__init__(x, y, 1, -1)

class CoinSuperieurDroit(Coin):
    def __init__(self, x, y):
        super().__init__(x, y, -1, -1)

## surfaces

# surface
# i1, i2 et i3 sont les index des blocs en contact avec le bloc élem
class Surface(Point):
    def __init__(self, x, y, i1, i2, i3):
        super().__init__(x, y)
        self.i1 = i1
        self.i2 = i2
        self.i3 = i3
    
    def getEquation(self):
        A = np.zeros(size**2)
        index = size * self.x + self.y
        A[index + self.i1] = k
        A[index + self.i2] = k
        A[index + self.i3] = k
        A[index] = -3 * k - h * dx
        return A
    
    def getConvectionTerm(self):
        return -h * dx * Tinf

class SurfaceVerticaleHaut(Surface):
    def __init__(self, x, y):
        super().__init__(x, y, -1, 1, -size)

class SurfaceVerticaleBas(Surface):
    def __init__(self, x, y):
        super().__init__(x, y, -1, 1, size)

class SurfaceHorizontaleDroite(Surface):
    def __init__(self, x, y):
        super().__init__(x, y, size, -size, 1)

class SurfaceHorizontaleGauche(Surface):
    def __init__(self, x, y):
        super().__init__(x, y, size, -size, -1)

       