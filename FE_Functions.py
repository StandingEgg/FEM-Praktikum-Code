import numpy as np
from Numerik_Tool import *
from matplotlib.pyplot import figure, plot,show, title


# Liefert die gewuenschten n Ansatzfunktionen
def getAnsatzfkt ():

    # lineare Ansatzfunktionen
    if N_num == 2:
        N_1 = lambda xi: 1./2 * (1-xi)
        N_2 = lambda xi: 1./2 * (1+xi)
        N = np.array ([N_1,N_2])                # abgelegt in einem Array

        # Nach xi abgeleitete Ansatzfunktionen
        # TODO Praktikumsaufgabe 4: - Ableitung der Ansatzfunktionen nach xi
        #                - Hilfe-Schlagwort: Mathematische Funktionen
        N_1_xi = -1./2
        N_2_xi =  1./2
        N_xi = np.array ([N_1_xi,N_2_xi])       # abgelegt in einem Array

    # quadratische Ansatzfunktionen
    if N_num == 3:
        # TODO Praktikumsaufgabe 5: - Definition der Ansatzfunktionen
        N_1
        N_2
        N_3
        N = np.array ([N_1,N_2,N_3])

        # Nach xi abgeleitete Ansatzfunktionen
        # TODO Praktikumsaufgabe 5: - Ableitung der Ansatzfunktionen nach xi
        N_1_xi
        N_2_xi
        N_3_xi
        N_xi = np.array ([N_1_xi,N_2_xi,N_3_xi])

    return N, N_xi

# Randbedingungen und Parameter
def defElements():


    # TODO Praktikumsaufgabe 4: - Definieren Sie hier Randbedingungen und Parameter wie
    #         z.B. Randkraft, Querschnittsflaeche, E-Modul...

    # NEUMANN-RB
    f_V = 0.1                               # Volumenlast
    F_R = 22500                               # Kraft am rechten Stabende

    # Geometrische-RB
    L = 1100                                 # Gesamtlaenge
    A = 225                                 # Querschnittsflaeche

    # Materialparameter
    E =2.1E5                                 # E-Modul


    #####  Werden erst in Termin 4 und 5 von Bedeutung sein => zunaechst mit Null initialisiert
    # TODO Praktikumsaufgabe 5: Materialparameter fuer isotrope Verfestigung definieren
    # TODO Praktikumsaufgabe 6: Materialparameter fuer kinematische Verfestigung definieren
    sigma_y   = 0                       # Fliessspannung
    K_elplast = 0                       # elastoplast. Tangentenmodul
    epsilon_p = 0                       # initiale plastische Dehnung
    m = 0                               # Exponent der isotropen Verfestigung
    H = 0                               # kinematische Verfestigung
    q = 0                               # initiale kinematische Verfestigungsvariable: q_0
    #####

    # Groesse des Arrays fuer plast Dehnung und kinem Verfestigung an die Anzahl von Gausspunkten anpassen
    if g_num == 1:
        epsilon_p_array = np.array([epsilon_p])
        q_array = np.array([q])

    elif g_num == 2:
        q_array = np.array([q,q])
        epsilon_p_array = np.array([epsilon_p,epsilon_p])
    else:
        print ("Mehr als 2 Gausspunkte werden nicht unterstuetzt!")

    # Definiere Elemente
    # TODO Praktikumsaufgabe 4: - Jedes Element kennt diverse Eigenschaften von sich selbst. Hier: globale Position des ersten Knotens, Laenge des Elements, ...
    #                - Vervollstaendigen Sie fehlende Eintraege fuer e2, e3 und e4
    #                - Hilfe-Schlagwort: Dictionary

    e1 = {'1.Knoten':0,       'h':2.*L/11, 'f_V':f_V, 'F_R':0,   'Position_F_R':0, 'E':E, 'A':A, 'sigma_y':sigma_y, 'K_elplast':K_elplast, 'm':m, 'H':H, 'q':q_array.copy(), 'epsilon_p':epsilon_p_array.copy() }
    e2 = {'1.Knoten':2.*L/11, 'h':2.*L/11, 'f_V':f_V, 'F_R':0,   'Position_F_R':0, 'E':E, 'A':A, 'sigma_y':sigma_y, 'K_elplast':K_elplast, 'm':m, 'H':H, 'q':q_array.copy(), 'epsilon_p':epsilon_p_array.copy() }
    e3 = {'1.Knoten':4.*L/11, 'h':2.*L/11, 'f_V':f_V, 'F_R':0,   'Position_F_R':0, 'E':E, 'A':A, 'sigma_y':sigma_y, 'K_elplast':K_elplast, 'm':m, 'H':H, 'q':q_array.copy(), 'epsilon_p':epsilon_p_array.copy() }
    e4 = {'1.Knoten':6.*L/11, 'h':1.*L/11, 'f_V':f_V, 'F_R':0,   'Position_F_R':0, 'E':E, 'A':A, 'sigma_y':sigma_y, 'K_elplast':K_elplast, 'm':m, 'H':H, 'q':q_array.copy(), 'epsilon_p':epsilon_p_array.copy() }
    e5 = {'1.Knoten':7.*L/11, 'h':4.*L/11, 'f_V':0,   'F_R':F_R, 'Position_F_R':1, 'E':E, 'A':A, 'sigma_y':sigma_y, 'K_elplast':K_elplast, 'm':m, 'H':H, 'q':q_array.copy(), 'epsilon_p':epsilon_p_array.copy() }
    return np.array([e1,e2,e3,e4,e5])

def Dirichlet_RB ():

    D = np.zeros(e_num*(N_num-1)+1)             # Verschiebungsfreiheitsgrade zunaechst mit Null belegt (willkuerlich).
    D[0] = 0                                    # An der ersten Stelle (Einspannung) zwingend Null.

    # TODO Praktikumsaufgabe 5: Verschiebungsrandbedinung am rechten Rand vorgeben

    return D


def Steifigkeitsmatrix(elemente, N_xi):

    # Ermitteln der Groesse der Steifigkeitsmatrix
    S_num = e_num*(N_num-1) + 1

    # Ermitteln der Steifigkeitsmatrix S
    S = np.zeros ((S_num,S_num))
    # TODO Praktikumsaufgabe 4:

    E_modul_list = []
    A_list = []
    e_length_list = []
    for i in range(e_num):
        e = defElements()
        E_modul = e[i]['E']
        E_modul_list.append(E_modul)
        A = e[i]['A']
        A_list.append(A)
        e_length = e[i]['h']
        e_length_list.append(e_length)


    for e in range(e_num):
        el = (N_num-1)*e
        N_1_xi = -1./2
        N_2_xi =  1./2
        N_xi = [N_1_xi, N_2_xi]
        for i, p in enumerate(N_xi):
            for j, q in enumerate(N_xi):
                Se = np.zeros((N_num, N_num))
                gauss_int = GaussQuad(lambda xi: p * q, g_num)
                Se_value = 2. / e_length_list[e] * E_modul_list[e] * A_list[e] * gauss_int
                Se[i][j] = Se_value

                for i in range(N_num):
                    for j in range(N_num):
                        S[el+i][el+j] += Se[i][j]

        # Schleife ueber die Elemente
            # Schleife ueber die Ansatzfkten der Testfunktion
                # Schleife ueber die Ansatzfkten der Formfunktion

                    # Berechnungsformel fuer die Steifigkeitsmatrix S. Die Assemblierung soll automatisch erfolgen.

                    # Hinweise zur Verwendung der Integrationsfunktion GaussQuad:
                    # Sollen zwei oder mehr Funktionen miteinander multipliziert und anschliessend integriert werden
                    # so wird dies folgendermassen uebergeben
                    # f_1 = lambda x: 2*x + 1
                    # f_2 = lambda x: x**2
                    # g_int = GaussQuad (lambda xi: f_1(xi)*f_2(xi),2)
    return S

def Lastvektor (elemente, N, N_xi, Lastschritt):

    # Ermitteln der Groesse des Lastvektors
    F_num = e_num*(N_num-1) + 1

    # Ermitteln des Kraftvektors F
    F = np.zeros ((F_num))

    A_list = []
    f_V_list = []
    e_length_list = []
    for i in range(e_num):
        e = defElements()
        A = e[i]['A']
        A_list.append(A)
        e_length = e[i]['h']
        e_length_list.append(e_length)
        f_V = e[i]['f_V']
        f_V_list.append(f_V)

    for e in range(e_num):
        el = (N_num-1)*e
        N_1 = lambda xi: 1./2 * (1-xi)
        N_2 = lambda xi: 1./2 * (1+xi)
        N = [N_1, N_2]

        for i, p in enumerate(N):
            fe = np.zeros((N_num))

            gauss_int = GaussQuad(lambda xi: p(xi) * f_V_list[i], g_num)
            fe_value = e_length_list[e] /2 * A_list[e] * gauss_int
            fe[i] = fe_value

            for i in range(N_num):

                F[el+i] += fe[i]



    # TODO Praktikumsaufgabe 4:
        # Schleife ueber die Elemente
            # Schleife ueber die Ansatzfkten der Testfunktion

                # Berechnungsformel fuer den Kraftvektor F. Die Assemblierung soll automatisch erfolgen.
    return F

def Reduzieren_Steifigkeit (S):

    # Groesse der Steifigkeitsmatrix
    S_num = S.shape[0]

    # Ermitteln des Kraftvektors unter Beruecksichtigung der Einspannung
    # TODO Praktikumsaufgabe 4:  - Steifigkeitsmatrix um die linke Einspannung reduzieren
    if Type_RB == 'NEUMANN-RB':
        S_red

    # Ermitteln des Kraftvektors wenn Verschiebungs-RB gesetzt
    # TODO Praktikumsaufgabe 5:  - Steifigkeitsmatrix um die linke und rechte Einspannung reduzieren
    elif Type_RB == 'DIRICHLET-RB':
        S_red

    return S_red

def Anpassen_Lastvektor (S, F, D):

    # Groesse des Lastvektors
    F_num = F.shape[0]

    # TODO Praktikumsaufgabe 4:  - Lastvektor um die linke Einspannung reduzieren
    if Type_RB == 'NEUMANN-RB':
        F

    # TODO Praktikumsaufgabe 5:  - Lastvektor um die linke Einspannung reduzieren
    #                 - Rechte Verschiebung mit Steifigkeitsmatrix verrechnen und Lastvektor damit anpassen
    elif Type_RB == 'DIRICHLET-RB':
        F

    return F

# Anpassen des Verschiebungsvektors
def Anpassen_Verschiebungsvektor(D, x):

    # Groesse des Verschiebungsvektors
    D_num = D.shape[0]

    # Passt den Verschiebungsvektor D um die berechneten Werte x an. Der erste Eintrag (linke Einspannung)
    # ist bereits enthalten.
    # TODO Praktikumsaufgabe 4:  - vervollstaendigen Sie die fehlende Codezeile
    if Type_RB == 'NEUMANN-RB':
        D

    # Passt den Verschiebungsvektor D um die berechneten Werte x an. Der erste (linke Einspannung)
    # und der letzte (Verschiebungs-RB) Eintrag sind bereits enthalten.
    # TODO Praktikumsaufgabe 5:  - vervollstaendigen Sie die fehlende Codezeile
    elif Type_RB == 'DIRICHLET-RB':
        D

    return D



# Elastischer Praediktor
def Elastischer_Praediktor (elemente, N_xi, D):

    # Anzahl der Gausspunkte
    if g_num == 1:
        g = np.array([0])
    elif g_num == 2:
        g = np.array([-(1./3)**0.5,+(1./3)**0.5])
    else:
        print ("Es sind nur 1-2 Gausspunkte moeglich!")

    # TODO Praktikumsaufgabe 4:  - Ziel ist es die elastische Materialantwort zu ermitteln
    #                  - Hierzu werden Dehnung und Spannung an den Gausspunkten ermittelt und in einem 1-2
    #                    dimensionalen Array abgespeichert
    #                  - Realisieren Sie eine Schleife ueber alle Elemente und berechnen sie die Spannungs-
    #                    und Dehnungswerte elementweise

    #Schleife ueber alle Elemente
    for e in range(e_num):

        # Dehnung (= epsilon) = N_xi*d;
        elemente[e]['Dehnung']

        # Spannung = E * (epsilon - epsilon^p)
        elemente[e]['Spannung']

    return

# Plastischer Korrektor
def Plastischer_Korrektor (elemente):

    if Type_Verf == 'ELASTISCH':
        # "plastisch" = False solange keine plastische Verfestigung auftritt
        plastisch = False

    elif Type_Verf == 'ISOTROP':
        # LINEARE ISOTROPE VERFESTIGUNG
        # Wir haben eine Verschiebungsrandbedingung (Dirichlet-RB); die abzutragende Dehnung ist vorgegeben.

        # "plastisch" = False solange keine plastische Verfestigung auftritt
        plastisch = False

        # TODO Praktikumsaufgabe 5:  - Ermitteln Sie die Fliessfunktion, ueberpruefen Sie ob die Fliessbedingung erfuellt ist oder nicht
        #                 - Ermitteln Sie die resultierende plastische Dehnung wenn die Fliessbedingung erfuellt ist
        #                 - Setzen Sie die Variable "plastisch" auf True sobald eine plastische Dehnung erfolgte

        # Schleife ueber die Elemente
            # Pruefen der Fliessbedingung
            # Schleife ueber die Gausspunkte
                    # Plastischer Korrektor


    elif Type_Verf == 'KINEMATISCH':
        # KINEMATISCHE VERFESTIGUNG
        # Wir haben eine Kraftrandbedingung (NEUMANN-RB); der Querschnitt des Dehnstabs bleibt konstant => die abzutragende Spannung ist vorgegeben.

        # "plastisch" = False solange keine plastische Verfestigung auftritt
        plastisch = False

        # TODO Praktikumsaufgabe 6:  - Ermitteln Sie die Fliessfunktion, ueberpruefen Sie ob die Fliessbedingung erfuellt ist oder nicht
        #                 - Ermitteln Sie die resultierende plastische Dehnung wenn die Fliessbedingung erfuellt ist
        #                 - Setzen Sie die Variable "plastisch" auf True sobald eine plastische Dehnung erfolgte

        # Schleife ueber die Elemente
            # Pruefen der Fliessbedingung
            # Schleife ueber die Gausspunkte
                    # Plastischer Korrektor


    return plastisch

def Postprocessing (elemente, S, F, D, N, N_xi):

    print ("\nSteifigkeitsmatrix S:\n"+str(S))
    print ("\nKraftvektor F:\n"+str(F))
    print ("\nVerschiebungsvektor D:\n"+str(D))

    # Dient der plot-funktion. Unterteilt den Bereich [-1,1] in 21 Teilbereiche
    xi = np.linspace(-1,1,21)

    # TODO Praktikumsaufgabe 4:  - Vervollstaendigen Sie die Funktion Postprocessing so, dass
    #                   Dehnungs- und Spannungsverlauf ueber die globale x Koordinate ausgeben werden.

    # Schleife ueber alle Elemente
    for e in range(e_num):
        elemente[e]['Verschiebungsfeld']    # Verschiebungsfeld = N*d
        elemente[e]['Dehnungsfeld']         # Dehnungsfeld (^= epsilon) = N_x*d;
        elemente[e]['epsilon_p_feld']    = make_func_of_epsilon_p(elemente[e]['epsilon_p'])(xi)*np.ones(21)
        elemente[e]['Spannungsfeld']        # Spannungsfeld = E * (epsilon - epsilon^p)


    # TODO Praktikumsaufgabe 4  - Elementweise ermittelte Felder (siehe oben) in einem Array zusammenfassen

    Spannungsfeld = []; Verschiebungsfeld = []; epsilon_p_feld = []; X = []
    # Schleife ueber alle Elemente
    for e in range(e_num):
        Spannungsfeld       # Element-Spannungsfelder aneinanderhaengen
        Verschiebungsfeld   # Element-Verschiebungsfelder aneinanderhaengen
        epsilon_p_feld      # Element-epsilon_p_felder aneinanderhaengen
        X                   # Transforamtion auf globale X-Koordinate


    # Ausgabe der Ergebisse in Form von Plots
    figure(1); title('Spannung', color='b', fontsize=20);          plot(X,Spannungsfeld)
    figure(2); title('Verschiebung', color='b', fontsize=20);      plot(X,Verschiebungsfeld)
    figure(3); title('plast. Dehnung', color='b', fontsize=20);    plot(X,epsilon_p_feld)
    show()

# Erstellt eine lineare mathematische Funktion aus dem Array der plastischen Dehnung
# Muss nicht mehr veraendert werden
def make_func_of_epsilon_p(element_epsilon_p):

    if element_epsilon_p.shape[0] == 1:
        f = lambda xi: element_epsilon_p

    elif element_epsilon_p.shape[0] == 2:
        f = lambda xi: (element_epsilon_p[1]-element_epsilon_p[0]) * 3.**0.5 / 2  * xi  +  (element_epsilon_p[1]+element_epsilon_p[0])/2.

    else:
        print ("Leider hat sich ein Fehler in der Funktion ""make_func_of_epsilon_p"" ereignet!")
        print ("Die Array-Groesse der plastischen Dehnung wird nicht unterstuetzt!")
        return 1
    return f

# Muss nicht mehr veraendert werden
def defglobalVariable (arg1, arg2, arg3, arg4, arg5):
    # arg1 = Anzahl an Gausspunkten
    # arg2 = Anzahl an Ansatzfunktionen
    # arg3 = Anzahl an Elementen
    # arg4 = Art der Verfestigung
    # arg5 = Art der rechten Randbedingung

    # definiert globale Variablen
    global g_num; global Type_Verf; global Type_RB; global N_num; global e_num

    g_num     = arg1
    N_num     = arg2
    e_num     = arg3
    Type_Verf = arg4
    Type_RB   = arg5