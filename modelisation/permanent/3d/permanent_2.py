import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from math import *
from matplotlib import cm
size = 8

# température avec échange par convection avec l'extérieur
Tinf = 300

# delta x correspondant à la longueur de controle
dx = .5

# conductivité thermique du matériau
k = 80

# puissance volumique
q = 1000

# coefficient de transfert par convection
h = 8


#Ici 0-> il n'y a rien
#1-> Il y a un bloc qui ne recoit pas de puissance volumique
#2-> Il y a un bloc qui recoit une puissance volumique


#C est la Matrice en 3D representant la structure avec des 0 si il y a rien et des 1 sinon
#Pour ce cas on prend tous les coeficients egaux à 1 sauf ceux sur les bords.
A=np.zeros((size,size))
for i in range(1,size-1) :
    for j in range(1,size-1) :
        A[i][j]+=1
for i in range(3,size-3) :
    for j in range(3,size-3) :   
        A[i][j]+=0
C=[A for k in range(size)]

#Les matrices A representent un coupe de l'espace 


#La convention choisie et que chaque matrice correspond à une coupe horizontale de l'espace
#On regarde l'objet du "dessus" et l'indice dans la liste qui continet les matrices correspond a la hauteur.



#Cette fonction sera utile par la suite pour trouver le minimum de la fonction
#Car on doit prendre au tout début une valeur non nulle du minimum(sinon min=0)
#Pas très important
def recherche_1er_coef_non_nul(A) :
    k=0 
    while A[k][0]==0 and k<(len(A)-1):
        k+=1
    return(k)
        
    

#on instancie le point a partir des coordonnées x,y,z
def getnature(x,y,z) :
    #On regarde si il y a des "blocs" autour cad si il y a des 1 ou des 2:
    c=[]
    #Si on se place aux extremités on a une erreur d'index c'est pourquoi
    #On vérifie avant qu'on est pas aux extremités a chaque fois.
    if x<(size-1) and (C[z][x+1][y]==1 or C[z][x+1][y]==2) :
       c+=[1]
    if x>0 and (C[z][x-1][y]==1 or C[z][x-1][y]==2):
       c+=[-1]
    if y<(size-1) and (C[z][x][y+1]==1 or C[z][x][y+1]==2):
       c+=[size]
    if y>0 and (C[z][x][y-1]==1 or C[z][x][y-1]==2):
       c+=[-size] 
    if z<(size-1) and (C[z+1][x][y]==1 or C[z+1][x][y]==2):
       c+=[size**2] 
    if z>0 and (C[z-1][x][y]==1 or C[z-1][x][y]==2) :
       c+=[-size**2]
    a=Point(x,y,z,c)
    return(a)                       
                    
def solve() :
    
    #La matrice A contient la matrice carrée utilisée pour le pivot de Gauss
    A=np.zeros((size**3,size**3))
    #La matrice B contient la matrice colonne utilisée pour le pivot
    B=np.zeros((size**3,1))
    
    for k in range(size**3) :
        #coordonnées du point
        z=k//(size**2)
        y=(k-z*size**2)//size
        x=k%size
        #Si il y a bien un"bloc" cet endroit on résout de manière normale
        if C[z][x][y]==1 :
            a=getnature(x,y,z)
            A[k]=a.getEquation()
            B[k]=a.getConvectionTerm()
        elif C[z][x][y]==2 :
            a=getnature(x,y,z)
            A[k]=a.getEquation()
            B[k]=a.getConvectionTerm()-q*(dx)**2
        #Si il n'y a rien a cet endroit on doit quand meme mettre au moins un coef non-nul
        #Sinon le pivot de Gauss bug on met donc le coef correspondant aux coordonées 
        #la temperature de ce point sera donc comptée comme nulle
        else :
            A[k][k]=1
            B[k]=0
    print(C)
    
    S=np.linalg.solve(A,B)
   
    return(S)
    
def visio_3D() :
    cmap=cm.coolwarm
    ax=Axes3D(plt.figure())
    S=solve()
    #On cherche a connaitre le min et le max pour bien exploiter le gradoient de couleur
    #Comme dit avant on prend le premier coef non nul pour le min sinon min=0
    j=recherche_1er_coef_non_nul(S)
    min=S[j][0]
    max=S[0][0]
    for k in range(len(S)) :
        if S[k][0] >max :
            max=S[k][0]
        if S[k][0]!=0 and S[k][0]<min :
            min=S[k][0]
    for x in range(size) :
        for y in range(size) :
            for z in range(size) :
                k=x+size*y+size**2*z
                #On trace une boule seulement si il y a un bloc
                if S[k][0]!=0 :
                    Theta,Phi=np.mgrid[0:pi:0.5,0:2*pi:1]
                    X=np.sin(Theta)*np.cos(Phi)*0.2+x
                    Y=np.sin(Theta)*np.sin(Phi)*0.2+y
                    Z=np.cos(Theta)*0.2+z
                    #Normallement 0-> bleu et 1-> foncé ici ont fait en sorte
                    #que les valeurs soient comprises dans cet intervalle
        
                    coef_couleur=(S[k][0]/(max-min)-min/(max-min))
                    ax.plot_surface(X,Y,Z,color=cmap(coef_couleur))          
    plt.title("Chaleur  maximale: "+str(int(max))+"K Chaleur minimale: "+str(int(min))+"K")
    plt.show()

def coupeHoriz(z) :
    ax=Axes3D(plt.figure())
    S=solve()
    T=np.zeros(size**2)
    #On prend seulement une matrice qui correspond aux valeurs voulues pour la coupe
    for k in range(size**2) :
        T[k]=S[z*size**2+k]
    X,Y =np.mgrid[0:size:1,0:size:1]
    #On fait en sorte  que les blocs ou il y a rien aient la temperature Tinf
    #Le graphique est moins ecrasé grace a ca
    Z=Tinf*np.ones((size,size))
    for i in range(size) :
        for j in range(size) :
            if T[i+j*size]!=0 :
                Z[i][j]=T[i+j*size]
    surf=ax.plot_surface(X,Y,Z,cmap=cm.coolwarm)
    
    plt.colorbar(surf)
    plt.show()











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
        return(A)
    # terme de convection à droite de l'égalité, il ne faut pas oublier
    # de lui sommer le terme de puissance volumique.
    def getConvectionTerm(self):
        return(- (6 - len(self.c)) * h * dx * Tinf)

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