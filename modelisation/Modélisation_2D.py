# -*- coding: utf-8 -*-
import numpy as np
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def modélisation_2D() :
    ax=Axes3D(plt.figure())
    #Voir la feuille sur le drive pour mieux comprendre le programme avec l'analytique
    dx=0.5 #coté d'un bloc en metre
    L=0.2 #profondeur du bloc
    h=10 #Coeficient de convection de l'air pour une paroi verticale (Cf Wikipedia)
    q=100 #Puissance volumique recue par le micro processeur
    G=80 #conductivité du matériau lambda (80 pour le fer)
    T=np.zeros((4,5)) #Matrice permettant de stocker les données relatives a la temperature
    T0=300 #température de l'air ambiant (32C)
    #Dans cet algorithme on se place dans le cas du régime permanent
    A=np.array([[-2*G-2*h*dx,G,0,G,0,0],[G,-3*G-h*dx,G,0,G,0],[0,G,0,G,-3*G,G],[G,0,0,-2*G-h*dx,G,0],[0,2*G,-3*G-h*dx,0,0,G],[0,0,G,0,2*G,-3*G]])
    B=np.array([[-q*(dx)**2-2*h*dx*T0],[-q*(dx)**2-h*dx*T0],[-q*(dx)**2],[-q*(dx)**2-h*dx*T0],[-q*(dx)**2-h*dx*T0],[-q*(dx)**2]])
    C=post_traitement(A,B)
    for j in range(3) :
        T[0][j]=C[j]
    for j in range(3,5) :
        T[0][j]=C[4-j]
    for j in range(3) :
        T[1][j]=C[j+3]
    for j in range(3,5) :
        T[1][j]=C[7-j]
    for j in range(3) :
        T[3][j]=C[j]
    for j in range(3,5) :
        T[3][j]=C[4-j]
    for j in range(3) :
        T[2][j]=C[j+3]
    for j in range(3,5) :
        T[2][j]=C[7-j]
    X,Y=np.mgrid[0:4:1,0:5:1]
    ax.plot_wireframe(X,Y,T)
    plt.show()
#Algorithme Pivot de Gauss qu'on doit utiliser par la suite
#Pivot de Gauss


def matrice_augmentée(A,X) :
    B=np.zeros((len(A),len(A[0])+1))
    for k in range(len(A)) :
        for j in range(len(A[0])) :
            B[k][j]=A[k][j]
        B[k][len(A[0])]=X[k]
    return(B)

def changement_ligne(A,i,j) :
    for k in range(len(A[i])) :
        (A[i][k],A[j][k])=(A[j][k],A[i][k])
    return(A)

def copie_profonde(A) :
    B=np.zeros((len(A),len(A[0])))
    for k in range(len(A)) :
        for l in range(len(A[0])) :
            B[k][l]=A[k][l]
    return(B)

def transvection(A,i,j,a) :
    for k in range(len(A[i])) :
          A[i][k]+=a*A[j][k]
    return(A)

def chercher_pivot(A,l) :
    b=0
    c=0
    for k in range(l,len(A)) :
        if abs(A[k][l])>b:
            b=abs(A[k][l])
            c=k
    return(c)

def dilatation(A,i,b) :
    for k in range(len(A[i])):
        A[i][k]=A[i][k]*b
    return(A)

def pivot_Gauss(A,X) :
    B=copie_profonde(A)
    B=matrice_augmentée(B,X)
    for k in range(len(B)) :
        i=chercher_pivot(B,k)
        B=changement_ligne(B,i,k)
        dilatation(B,k,1/B[k][k])
        for c in range(k+1,len(B)) :
            B=transvection(B,c,k,-B[c][k])   #On met des 0 partout sur la colonne
        print(B)
    return(B)

def max(A,a) :
    m=0
    for k in range(a) :
        for j in range(len(A[k])) :
            if A[k][j]>m :
                m=A[k][j]
    return(m)

def post_traitement(A,X) :
    if len(A)!=len(X) :
        return("dimension incompatible")

    B=pivot_Gauss(A,X)
    a=len(B)
    b=len(B[0])
    for k in range(a) :
        l=a-k-1
        for i in range(l) :
            B=transvection(B,i,l,-B[i][l])
    print(B)
    if max(B,a)>(10**10)*max(A,a) :
        return("resultat incohérent")
    Sol=np.zeros((len(B),1))
    for k in range(len(B)) :
        Sol[k]=B[k][b-1]
    return(Sol)
