# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pyplot as plt
import copy

# ------------- Paramètres physiques -------------------------------------------------------------------

a = 0.06 # Dimension de l'arrete du cubre de platre(m)
l = 0.01 # Largeur de la cavité (m)(dimension selon x)
L = 0.01 # Longueur de la cavité (m) (dimension selon z)
lambd = 0.6 #lambda = 0.6 # Conductivité thermique (W/m/K)
h = 50 # Coefficient de convection (W/m2/K)
Tint = 35# Température à l'intérieur de la cavité (°C)
Text = 20 # Température ambiante (°C)

# ------------- Paramètres numériques ------------------------------------------------------------------

n = 13 # Nombre de pas de discrétisations (nombre pair)
epsilon = 1e-4 # Critère d'arret
d = a/n # Dimension du pas de discrétisation



# ------------- Définition des indices limites de la cavité : fonction CAVITE ---------------------------

def Cavite (a, l, L, d):
    i_min = int((a/(2*d))-(l/d))+1
    i_max = int((a/(2*d))+(l/d))
    j_min = int((a/(2*d))-(L/d))+1
    j_max = int((a/(2*d))+(L/d))
    A = np.array([[float(Text) for i in range(n+1)] for j in range(n+1)])
    for i in range(i_min, i_max+1):
        for j in range(j_min, j_max+1):
            A[i,j]=float(Tint)
    return A

# ------------- Calcul des températures par méthode itérative : fonction JACOBI -------------------------
def jacobis(T0,n,epsilon,Tint):
    incr =100
    T=np.array(T0,copy=True)
    while (incr>epsilon):
        T_old =copy.deepcopy(T)
        for i in range (1,n):
            for j in range (1,n):
                if not T[i,j]==Tint:
                    T[i,j]=1/4*(T[i+1,j]+T[i-1,j]+T[i,j-1]+T[i,j+1])
        incr = np.linalg.norm(T-T_old)
        for i in range (0,n+1):
            T[n,i]=T[n-1,i]
            T[i,n]= T[i,n-1]
            T[i,0]=T[i,1]
            T[0,i]=T[1,i]/(1+h*d/lambd)

    return T


# ------------- Calcul du flux par Simpson 1/3 : fonction FLUX2D -----------------------------------------



# ------------------------ PROGRAMME PRINCIPAL ----------------------------------------------------------
T = np.zeros((n+1,n+1),float)
T0=Cavite(a, l, L, d)

T=jacobis(T0,n,epsilon,Tint)
