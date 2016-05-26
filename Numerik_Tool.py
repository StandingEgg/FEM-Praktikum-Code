import numpy as np
import math
from FE_Functions import *





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
        results = w * f(p)


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
