# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:33:14 2016

@author: Yifei
"""

import numpy as np
import math
from FE_Functions31 import *

# Gausssche Quadratur
def GaussQuad (f, g_num):
    # f=zu integrierende Funktion, g_num=Anzahl der Stuetzpunkte, liefert den exakten Integrationswert fuer ein Polynom vom Grad (2*g_num-1)

    # TODO Praktikumsaufgabe 1: - Ueberpruefen Sie welche Anzahl an erforderlichen Stuetzpunkten uebergeben wurde
    #                - Definieren Sie sich fuer die jeweilige Anzahl an Stuetzpunkten
    #                  die Position der Stuetzpunkte und die Gewichtung
    #                - Werten Sie die Summe c_i * f(x_i) aus und liefern Sie das Ergebnis mit dem "return" Befehl zurueck


    #dX_dxi = (b-a)/2.0
    if g_num == 1:
        points = [0]
        weights = [2]
    elif g_num == 2:
        points = [-np.sqrt(1/3.0),np.sqrt(1/3.0)]
        weights = [1,1]
    elif g_num == 3:
        points = [-np.sqrt(3/5.0),0,np.sqrt(3/5.0)]
        weights = [5/9.0,8/9.0,5/9.0]

    results = 0
    for i, p in enumerate(points):
        w = weights[i]
        #x = (p+1)/2.0*(b-a)+a
        #results += w * f(x) * dX_dxi
        # Calculate with f_xi, not f_x !
        results += w * f(p)


    return results


f_x = lambda x: math.sin(2*x)

f_xi = lambda xi: math.sin(2*((xi+1)/2.0*(b-a)+a)*dX_dxi)


# Vorwarts-Rueckwaerts-Einsetzen
def Vorwaerts_Rueckwaerts_Einsetzen(L,R,b):
    y = np.linalg.solve(L, b)
    x = np.linalg.solve(R, y)

    # Vorwaertseinsetzen:  # Ly = b  => y
    # TODO Praktikumsaufgabe 4: - Loesen Sie das Gleichungssystem durch Vorwaertseinsetzen
    #                           - Zur Kontrolle: Ly = b  => y = np.linalg.solve(L,b)

    # Rueckwaertseinsetzen:  Rx = y  => x
    # TODO Praktikumsaufgabe 4: - Loesen Sie das Gleichungssystem durch Rueckwaertseinsetzen
    #                           - Zur Kontrolle: Rx = y  => x = np.linalg.solve(R,y)
    return x


# LR-Zerlegung
# Muss nicht mehr veraendert werden
def LR_Zerlegung (S):
    # Angabe: Ax = b
    # Ziel:   Ly = b  => y
    #         Rx = y  => x
    #================= Ermittlung der oberen und unteren Dreiecksmatrix

    A = S.copy()
    n = A.shape[0]
    L = np.zeros((n, n))
    R = np.zeros((n, n))
    for i in range(n):
        # Bestimmen von R
        for j in range(i,n):
            for k in range(i):
                A[i][j] -= A[i][k] * A[k][j]

        # Bestimmen von L
        for j in range (i+1, n):
            for k in range (i):
                A[j][i] -= A[j][k] * A[k][i]
            A[j][i] /= A[i][i]

    for i in range(n):
        for j in range (i,n):
            R[i][j] = A[i][j]

    L=A-R
    for i in range(n):
        L[i][i]=1

    return L,R

# def Newton_Iteration
def dx(f, x):
    return abs(0-f(x))

def newtons_method(f, df, x0, e):
    delta = dx(f, x0)
    while delta.all() > e:
        x0 = x0 - f(x0)/df(x0)
        delta = dx(f, x0)
    return x0

    # print ('Root is at: ', x0)
    # print ('f(x) at root is: ', f(x0))

def derivative(f):
# Computes the numerical derivative of a function.
    def df(x, h=0.1e-5):
        return ( f(x+h/2) - f(x-h/2) )/h
    return df
