import numpy as np
import math

a = math.pi/2.0
b = math.pi
dX_dxi = (b-a)/2.0
# Gausssche Quadratur
def GaussQuad (f, g_num):
    # f=zu integrierende Funktion, g_num=Anzahl der Stuetzpunkte, liefert den exakten Integrationswert fuer ein Polynom vom Grad (2*g_num-1)

    # TODO Praktikumsaufgabe 1: - Ueberpruefen Sie welche Anzahl an erforderlichen Stuetzpunkten uebergeben wurde
    #                - Definieren Sie sich fuer die jeweilige Anzahl an Stuetzpunkten
    #                  die Position der Stuetzpunkte und die Gewichtung
    #                - Werten Sie die Summe c_i * f(x_i) aus und liefern Sie das Ergebnis mit dem "return" Befehl zurueck



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
        x = (p+1)/2.0*(b-a)+a
        results += w * f(x) * dX_dxi


    return results


f_x = lambda x: math.sin(2*x)

f_xi = lambda xi: math.sin(2*((xi+1)/2.0*(b-a)+a)*dX_dxi)

F_analytisch = -1.0
F_ein_stuetzpunkt = GaussQuad(f_x, 1)
F1_relative_fehler = abs((F_ein_stuetzpunkt - F_analytisch)/F_analytisch)
# relative fehler bei einem Stuetzpunkt: 0.5707963267948966
F_zwei_stuetzpunkt = GaussQuad(f_x, 2)
F2_relative_fehler = abs((F_zwei_stuetzpunkt - F_analytisch)/F_analytisch)
F_drei_stuetzpunkt = GaussQuad(f_x, 3)
# relative fehler bei zwei Stuetzpunkten: 0.03209021267443135
F3_relative_fehler = abs((F_drei_stuetzpunkt - F_analytisch)/F_analytisch)
# relative fehler bei drei Stuetzpunkten: 0.0006944568038718124
