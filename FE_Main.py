import numpy as np
from FE_Functions import *

# Geometrie
#            f
#  -> -> -> -> -> -> -> ->
# 0==x==0==x==0==x==0=x=0====x====0 ---> F_R
#  -->(1)=0                                   globale Verschiebungsfreiheitsgrade
#        -->(2)
#              -->(3)
#                   -->(4)
#                        -->(5)
#                                 -->(6)

#=============================================  PREPROCESSING  =========================================================

# -------------------------------
# Diese Variablen dienen als Steuerelemente und koennen/muessen geaendert werden
# Anzahl an Gausspunkten
g_num = 2

# Anzahl an Ansatzfunktionen
N_num = 3

# Art des Materialverhaltens.       ELASTISCH, ISOTROP oder KINEMATISCH
Type_Verf = 'KINEMATISCH'

# Art der rechten Randbedingung.    DIRICHLET-RB  oder NEUMANN-RB
Type_RB = 'NEUMANN-RB'
# ---------------------------------

# Anzahl an Elementen
e_num = 5

# Definiert globale Variablen fuer FE_Functions
defglobalVariable(g_num, N_num, e_num, Type_Verf, Type_RB)

# Ansatzfunktionen
N, N_xi = getAnsatzfkt ()

# Definiert Elemente
elemente = defElements()

# Erzeugt entsprechend den gesetzten DIRICHLET-RB einen Verschiebungsvektor
D = Dirichlet_RB ()

# Lastschritte
if Type_Verf == 'ISOTROP' or Type_Verf == 'ELASTISCH':
    # Das Lastschrittfeld besteht lediglich aus einem Eintrag => die Schleife weiter unten wird nur einmal mit dem
    # Wert 1 Durchlaufen
    Lastschrittfeld = np.array([1])

elif Type_Verf == 'KINEMATISCH':
    # -----------------------------------------------------------------------------------
    # Wird erst in Praktikumsaufgabe 5 von Interesse sein
    # Hier wird das Lastschrittfeld aufgespannt. Spaeter wird jeder Eintrag im Lastschrittfeld durchlaufen.
    Lastschrittfeld = np.arange(0,1.0,0.1); Lastschrittfeld = np.append(Lastschrittfeld,np.arange(+1,-1.0,-0.1)); Lastschrittfeld = np.append(Lastschrittfeld,np.arange(-1.0,+1.1,+0.1))

    # hier werden zu jedem Lastschritt die entspr. Werte fuers Postpro. abgespeichert.
    Dehnung = []; Spannung = []; kinVerfVar = []
    # -----------------------------------------------------------------------------------
else:
    print ("Fuer das Materialverhalten sind nur die Keywoerter ELASTISCH, ISOTROP und KINEMATISCH zulaessig!")

#=============================================  SOLVER  =============================================================

# Berechnung der Steifigkeitsmatrix
S = Steifigkeitsmatrix(elemente, N_xi)

# Ermitteln der Steifigkeitsmatrix und des Kraftvektors unter Beruecksichtigung der Einspannung
S_red = Reduzieren_Steifigkeit (S)

# Loesung des Gleichgewichts mit Hilfe der LR-Zerlegung und anschliessendem Vorwaerts-Rueckwaerts-Einsetzen
L,R = LR_Zerlegung(S_red)


# laeuft die Lastschritte ab. Bei isotropem oder elastischem Materialverhalten gibt es nur einen Lastschritt = 1
for Lastschritt in Lastschrittfeld:

    # Fuer den elastischen Fall wird die while-Schleife nur einmal durchlaufen!!
    # Im isotropen Fall wird die while-Schleife solange durchlaufen bis fuer die Fliessfunktion == 0 gilt
    # Gesteuert wird die while-Schleife mit Hilfe der boolschen Variable "plastisch"
    plastisch = True
    while plastisch:

        # Berechnung der Steifigkeitsmatrix
        F = Lastvektor(elemente,N, N_xi, Lastschritt)

        # Passt den Lastvektor nach der gewaehlten rechten Randbedingung an
        F = Anpassen_Lastvektor (S, F, D)

        # Loesung des Gleichgewichts
        x = Vorwaerts_Rueckwaerts_Einsetzen(L,R,F)

        # Erweitert den Verschiebungsvektor um die berechneten Werte. Der erste Eintrag (Einspannung) bleibt unberuehrt.
        D = Anpassen_Verschiebungsvektor(D, x)

        # Elastischer Praediktor
        Elastischer_Praediktor (elemente, N_xi, D)

        # -----------------------------------------------------------------------------------
        # Wird erst in Praktikumsaufgabe 4 und 5 von Interesse sein
        # Plastischer Korrektor
        plastisch = Plastischer_Korrektor (elemente)
        # -----------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------
    # Wird erst in Praktikumsaufgabe 5 von Interesse sein
    if Type_Verf == 'KINEMATISCH':
        Spannung    = np.append(Spannung,elemente[0]['Spannung'][0])     # zu jedem Lastschritt wird der Spannungswert an Element (0) abgespeichert
        Dehnung     = np.append(Dehnung,elemente[0]['Dehnung'][0])       # zu jedem Lastschritt wird der Dehnungswert an Element (0) abgespeichert
        kinVerfVar  = np.append(kinVerfVar,elemente[0]['q'][0])          # zu jedem Lastschritt wird die kin. Verf. Variable an Element (0) abgespeichert
    # -----------------------------------------------------------------------------------

##============================================  POSTPROCESSING  =========================================================




if Type_Verf == 'ISOTROP' or Type_Verf == 'ELASTISCH':
    Postprocessing(elemente, S, F, D, N, N_xi)

# -----------------------------------------------------------------------------------
# Wird erst in Praktikumsaufgabe 5 von Interesse sein
elif Type_Verf == 'KINEMATISCH':
    title('Spannungs-Dehnungs-Diagramm', color='b', fontsize=20)
    plot(Dehnung,Spannung)
    plot(Dehnung,kinVerfVar)
    show()
# -----------------------------------------------------------------------------------
