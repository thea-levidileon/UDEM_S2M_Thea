import numpy as np
import numpy.matlib
import sys
#import ipopt
# sys.path.append('/home/user/anaconda3/envs/qpoases/lib')
# sys.path.append('/home/user/anaconda3/envs/qpoases/include/casadi/build/lib')
import casadi as cas
# from casadi import MX, SX, sqrt
# import biorbd
import matplotlib.pyplot as plt
# from os.path import dirname, join as pjoin
# import scipy.io as sio
from IPython import embed
from mpl_toolkits import mplot3d

n=15 #nombre de mailles sur le grand cote
m=9#nombre de mailles sur le petit cote
Masse_centre=150
PLT_FLAG=0

Nb_ressorts=2*n*m+n+m #nombre de ressorts non obliques total dans le modele
Nb_ressorts_cadre=2*n+2*m #nombre de ressorts entre le cadre et la toile
Nb_ressorts_croix=2*(m-1)*(n-1) #nombre de ressorts obliques dans la toile
Nb_ressorts_horz=n * (m - 1) #nombre de ressorts horizontaux dans la toile (pas dans le cadre)
Nb_ressorts_vert=m * (n - 1) #nombre de ressorts verticaux dans la toile (pas dans le cadre)

L_toile = 2.134
L_cadre=2.134+0.35
dL = 2 * L_toile / (n - 1)
#Largeur
l_toile = 1.07
l_cadre=1.07+0.38
dl = 2 * l_toile / (m - 1)

ind_milieu=int((m * n - 1) / 2)
ind_milieu_droite = int((m * n - 1) / 2 - n)
ind_milieu_haut =  int((m * n - 1) / 2 + 1)
ind_milieu_gauche = int((m * n - 1) / 2 + n)
ind_milieu_bas =  int((m * n - 1) / 2 - 1)
ind_milieu_bas_droite=int((m * n - 1) / 2 - n - 1)
ind_milieu_haut_droite = int((m * n - 1) / 2 - n + 1)
ind_milieu_haut_gauche= int((m * n - 1) / 2 + n + 1)
ind_milieu_bas_gauche= int((m * n - 1) / 2 + n - 1)


def Points_ancrage () :
    # ancrage :
    Pt_ancrage = np.zeros((3, Nb_ressorts_cadre))
    # cote droit :
    for i in range(n):
        Pt_ancrage[:, i] = np.array([l_cadre, -L_toile + i * dL, 0])
    # cote haut :
    for j in range(n, n + m):
        Pt_ancrage[:, j] = np.array([l_toile - (j - n) * dl, L_cadre, 0])
    # cote gauche :
    for k in range(n + m, 2 * n + m):
        Pt_ancrage[:, k] = np.array([-l_cadre, L_toile - (k - m - n) * dL, 0])
    # cote bas :
    for h in range(2 * n + m, Nb_ressorts_cadre):
        Pt_ancrage[:, h] = np.array([-l_toile + (h - 2 * n - m) * dl, -L_cadre, 0])
    return Pt_ancrage

def Spring_bouts(Pt,Pt_ancrage):
    # Definition des ressorts (position, taille)
    Spring_bout_1 = []

    # RESSORTS ENTRE LE CADRE ET LA TOILE
    for i in range(0, Nb_ressorts_cadre):
        Spring_bout_1 = cas.horzcat(Spring_bout_1, Pt_ancrage[:, i])

    # RESSORTS HORIZONTAUX : il y en a n*(m-1)
    for i in range(Nb_ressorts_horz):
        Spring_bout_1 = cas.horzcat(Spring_bout_1, Pt[:, i])

    # RESSORTS VERTICAUX : il y en a m*(n-1)
    for i in range(n - 1):
        for j in range(m):
            Spring_bout_1 = cas.horzcat(Spring_bout_1, Pt[:, i + n * j])

    Spring_bout_2 = []

    # RESSORTS ENTRE LE CADRE ET LA TOILE
    for i in range(0, n):
        Spring_bout_2 = cas.horzcat(Spring_bout_2, Pt[:, i])  # points droite du bord de la toile
    for i in range(n - 1, m * n, n):
        Spring_bout_2 = cas.horzcat(Spring_bout_2, Pt[:, i])  # points hauts du bord de la toile
    for i in range(m*n-1,n * (m - 1)-1, -1):
        Spring_bout_2 = cas.horzcat(Spring_bout_2, Pt[:, i])  # points gauche du bord de la toile
    for i in range(n * (m - 1), -1, -n):
        Spring_bout_2 = cas.horzcat(Spring_bout_2, Pt[:, i])  # points bas du bord de la toile

    # RESSORTS HORIZONTAUX : il y en a n*(m-1)
    for i in range(n, n * m):
        Spring_bout_2 = cas.horzcat(Spring_bout_2, Pt[:, i])

    # RESSORTS VERTICAUX : il y en a m*(n-1)
    for i in range(1, n):
        for j in range(m):
            Spring_bout_2 = cas.horzcat(Spring_bout_2, Pt[:, i + n * j])

    return (Spring_bout_1,Spring_bout_2)

def Spring_bouts_croix(Pt):
    #RESSORTS OBLIQUES : il n'y en a pas entre le cadre et la toile

    Spring_bout_croix_1=[]
    #Pour spring_bout_1 on prend uniquement les points de droite des ressorts obliques
    for i in range ((m-1)*n):
        Spring_bout_croix_1=cas.horzcat(Spring_bout_croix_1,Pt[:,i])
        #a part le premier et le dernier de chaque colonne, chaque point est relie a deux ressorts obliques
        if (i+1)%n!=0 and i%n!=0 :
            Spring_bout_croix_1 = cas.horzcat(Spring_bout_croix_1, Pt[:, i])

    Spring_bout_croix_2=[]
    #Pour spring_bout_2 on prend uniquement les points de gauche des ressorts obliques
    #pour chaue carre on commence par le point en haut a gauche, puis en bas a gauche
    #cetait un peu complique mais ca marche, faut pas le changer
    j=1
    while j<m:
        for i in range (j*n,(j+1)*n-2,2):
            Spring_bout_croix_2 = cas.horzcat(Spring_bout_croix_2, Pt[:, i + 1])
            Spring_bout_croix_2 = cas.horzcat(Spring_bout_croix_2, Pt[:, i])
            Spring_bout_croix_2 = cas.horzcat(Spring_bout_croix_2, Pt[:, i+ 2])
            Spring_bout_croix_2 = cas.horzcat(Spring_bout_croix_2, Pt[:, i + 1])
        j+=1
    return Spring_bout_croix_1,Spring_bout_croix_2

def Energie_func(Pt, Pt_ancrage, k, M, l_repos,k_croix,l_repos_croix):

    Spring_bout_1, Spring_bout_2=Spring_bouts(Pt,Pt_ancrage)
    Spring_bout_croix_1, Spring_bout_croix_2=Spring_bouts_croix(Pt)

    Energie = cas.MX.zeros(1)
    for i in range(Nb_ressorts):
        Energie += 0.5 * k[i] * (cas.norm_fro(Spring_bout_2[:, i] - Spring_bout_1[:, i]) - l_repos[i])**2
    for i_croix in range (Nb_ressorts_croix):
        Energie += 0.5 * k_croix[i_croix] * (cas.norm_fro(Spring_bout_croix_2[:, i_croix] - Spring_bout_croix_1[:, i_croix]) - l_repos_croix[i_croix]) ** 2

    Energie += cas.sum2(M.T * 9.81 * Pt[2, :])

    func = cas.Function('Energie_func', [Pt], [Energie]).expand()

    return func

def ForceEquilib_func(Pt, Pt_ancrage, k, M, l_repos,k_croix,l_repos_croix):

    Spring_bout_1, Spring_bout_2 = Spring_bouts(Pt, Pt_ancrage)
    Spring_bout_croix_1, Spring_bout_croix_2 = Spring_bouts_croix(Pt)

    Vect_unit_dir_F = (Spring_bout_2 - Spring_bout_1) / cas.norm_fro(Spring_bout_2 - Spring_bout_1)
    F_spring = cas.MX.zeros(3, Nb_ressorts)
    for ispring in range(Nb_ressorts):
        F_spring[:, ispring] = Vect_unit_dir_F[:, ispring] * k[ispring] * (
                cas.norm_fro(Spring_bout_2[:, ispring] - Spring_bout_1[:, ispring]) - l_repos[ispring])

    Vect_unit_dir_F_croix = (Spring_bout_croix_2 - Spring_bout_croix_1) /cas.norm_fro(Spring_bout_croix_2 - Spring_bout_croix_1)
    F_spring_croix = cas.MX.zeros(3, Nb_ressorts_croix)
    for ispring in range(Nb_ressorts_croix):
        F_spring_croix[:, ispring] = Vect_unit_dir_F_croix[:, ispring] * k_croix[ispring] * (
                cas.norm_fro(Spring_bout_croix_2[:, ispring] - Spring_bout_croix_1[:, ispring]) - l_repos_croix[ispring])

    F_spring=cas.horzcat(F_spring,F_spring_croix)
    F_spring=F_spring_croix

    func = cas.Function('ForceEquilib_func', [Pt], [F_spring])  #.expand()

    return func

def ForceEquilib_centre_func(Pt, k, M, l_repos, Masse_centre,k_croix_tab,l_repos_croix):
    ind_milieu=int((m * n - 1) / 2)

    ind_milieu_droite = int((m * n - 1) / 2 - n)
    ind_milieu_haut =  int((m * n - 1) / 2 + 1)
    ind_milieu_gauche = int((m * n - 1) / 2 + n)
    ind_milieu_bas =  int((m * n - 1) / 2 - 1)

    ind_milieu_bas_droite=int((m * n - 1) / 2 - n - 1)
    ind_milieu_haut_droite = int((m * n - 1) / 2 - n + 1)
    ind_milieu_haut_gauche= int((m * n - 1) / 2 + n + 1)
    ind_milieu_bas_gauche= int((m * n - 1) / 2 + n - 1)


    k_4 = k[np.array([Nb_ressorts_cadre+ind_milieu_droite,
            Nb_ressorts_cadre + Nb_ressorts_horz + ind_milieu,
            Nb_ressorts_cadre+ ind_milieu,
            Nb_ressorts_cadre +  Nb_ressorts_horz + ind_milieu_bas])]
    l_repos_4 = l_repos[np.array([Nb_ressorts_cadre+ind_milieu_droite,
                        Nb_ressorts_cadre + Nb_ressorts_horz + ind_milieu,
                        Nb_ressorts_cadre+ ind_milieu,
                        Nb_ressorts_cadre +  Nb_ressorts_horz + ind_milieu_bas])]

    k_4_croix=k_croix_tab[np.array([ind_milieu_droite,
                          ind_milieu_haut,
                          ind_milieu_haut_gauche,
                          ind_milieu_bas_gauche])]
    l_repos_croix_4=l_repos_croix[np.array([ind_milieu_droite,
                                  ind_milieu_haut,
                                  ind_milieu_haut_gauche,
                                  ind_milieu_bas_gauche])]

    Spring_bout_1 = cas.horzcat( Pt[:, ind_milieu_droite], Pt[:, ind_milieu_haut ], Pt[:, ind_milieu_gauche], Pt[:, ind_milieu_bas])
    Spring_bout_2 = cas.horzcat(Pt[:, ind_milieu], Pt[:, ind_milieu],Pt[:, ind_milieu], Pt[:, ind_milieu])

    Spring_bout_croix_1=cas.horzcat( Pt[:, ind_milieu_bas_droite], Pt[:, ind_milieu_haut_droite], Pt[:, ind_milieu], Pt[:, ind_milieu])
    Spring_bout_croix_2=cas.horzcat( Pt[:, ind_milieu], Pt[:, ind_milieu],Pt[:,ind_milieu_bas_gauche],Pt[:,ind_milieu_haut_gauche])


    Vect_unit_dir_F = (Spring_bout_1 - Spring_bout_2) / cas.norm_fro(Spring_bout_1 - Spring_bout_2)
    F_spring = cas.MX.zeros(3, 4)
    for ispring in range(4):
        F_spring[:, ispring] = Vect_unit_dir_F[:, ispring] * k_4[ispring] * (
                    cas.norm_fro(Spring_bout_1[:, ispring] - Spring_bout_2[:, ispring]) - l_repos_4[ispring])

    Vect_unit_dir_F_croix = (Spring_bout_croix_2 - Spring_bout_croix_1) / cas.norm_fro(
        Spring_bout_croix_2 - Spring_bout_croix_1)
    F_spring_croix = cas.MX.zeros(3, 4)
    for ispring in range(4):
        F_spring_croix[:, ispring] = Vect_unit_dir_F_croix[:, ispring] * k_4_croix[ispring] * (
                cas.norm_fro(Spring_bout_croix_2[:, ispring] - Spring_bout_croix_1[:, ispring]) - l_repos_croix_4[ispring])

    Force_Masse = cas.MX.zeros(3)
    Force_Masse[2] = - (Masse_centre + M[ind_milieu]) * 9.81

    Force_tot = cas.sum2(F_spring) + cas.sum2(F_spring_croix) + Force_Masse
    # Force_tot = cas.sum2(F_spring_croix) + Force_Masse
    func = cas.Function('ForceEquilib_func', [Pt], [Force_tot])  #.expand()

    return func


def Optimisation_toile(Masse_centre, Pt, Pt_ancrage, k, M, l_repos, Pos_repos,k_croix_tab,l_repos_croix):
    Pt = cas.MX.sym('Pt', 3, n*m) #Pt est une variable que lon va utiliser dans les fonctions
    # Pt=cas.MX(Pt)
    Pt_ancrage = cas.MX(Pt_ancrage)
    k = cas.MX(k)
    M = cas.MX(M)
    l_repos = cas.MX(l_repos)
    k_croix_tab=cas.MX(k_croix_tab)
    l_repos_croix = cas.MX(l_repos_croix)

    #fonctions utilisant come arguments toutes les variables symboliques
    Energie = Energie_func(Pt, Pt_ancrage, k, M, l_repos,k_croix_tab,l_repos_croix)
    ForceEquilib_centre = ForceEquilib_centre_func(Pt, k, M, l_repos, Masse_centre,k_croix_tab,l_repos_croix)
    ForceEquilib = ForceEquilib_func(Pt, Pt_ancrage, k, M, l_repos,k_croix_tab,l_repos_croix)

    w = [] #vecteur de variables
    w0 = [] #conditions initiales
    lbw = []
    ubw = []
    g = [] #contraintes
    lbg = []
    ubg = []
    Pt_post = cas.MX.zeros(3, n*m)
    for i in range(n*m):
        Pt_var = cas.MX.sym(f'Pt_{i}', 3) #pour chaque point, Pt_var = fonction de Pt[i]
        Pt_post[:, i] = Pt_var
        w += [Pt_var]
        w0 += [Pos_repos[:, i]]
        lbw += [Pos_repos[0, i] - 3, Pos_repos[1, i] - 3, Pos_repos[2, i] - 3] # [np.zeros(3) - 0.00001]  #
        ubw += [Pos_repos[0, i] + 3, Pos_repos[1, i] + 3, Pos_repos[2, i] + 3] # [np.zeros(3)+ 0.00001]  #

    g += [ForceEquilib_centre(Pt_post)]
    lbg += [-10*np.ones(3)] # np.zeros(3) - 1e-10]
    ubg += [10*np.ones(3)] # np.zeros(3) + 1e-10]

    F_masses = cas.MX.zeros(3, n*m)
    F_masses[2, :] = - M * 9.81

    F_spring = ForceEquilib(Pt_post)

    obj = Energie(Pt_post) # + Eq_Obj

    qp = {'x': cas.vertcat(*w), 'f': obj, 'g':  cas.vertcat(*g)}
    solver = cas.nlpsol('solver', 'ipopt', qp) #chercher sur intenret doc casadi, chercher exeples ipopt
    sol = solver(x0=cas.vertcat(*w0), lbg=cas.vertcat(*lbg), ubg=cas.vertcat(*ubg), lbx=cas.vertcat(*lbw), ubx=cas.vertcat(*ubw))

    Solution = sol['x']

    Pt = np.zeros((3, n*m)) #la on fait sortir le vrai Pt, pas juste symbolique
    for j in range(n*m):
        Pt[:, j] = np.reshape(Solution[3*j : 3*j + 3], (3)) #Pt prend la valeur de la solution optimale

    return Pt


def Force_calc(Pt,Masse_centre):

    #PARAMETRES ET POINTS
    #Longueur
    L_toile = 2.134
    L_cadre=2.134+0.35
    dL = 2 * L_toile / (n - 1)
    #Largeur
    l_toile = 1.07
    l_cadre=1.07+0.38
    dl = 2 * l_toile / (m - 1)

    # ancrage :
    Pt_ancrage = np.zeros((3, Nb_ressorts_cadre))
    # cote droit :
    for i in range(n):
        Pt_ancrage[:, i] = np.array([l_cadre, -L_toile + i * dL, 0])
    # cote haut :
    for j in range(n, n + m):
        Pt_ancrage[:, j] = np.array([l_toile - (j - n) * dl, L_cadre, 0])
    # cote gauche :
    for k in range(n + m, 2 * n + m):
        Pt_ancrage[:, k] = np.array([-l_cadre, L_toile - (k - m - n) * dL, 0])
    # cote bas :
    for h in range(2 * n + m, Nb_ressorts_cadre):
        Pt_ancrage[:, h] = np.array([-l_toile + (h - 2 * n - m) * dl, -L_cadre, 0])

    # repos :
    Pos_repos = np.zeros((3, n*m))
    for i in range(n*m):
        Pos_repos[:, i] = np.array([l_toile - (i // n) * dl, -L_toile + (i % n) * dL, 0])

    # k trouves a partir du programme 5x3:
    k1 = (5 / n) * 3266.68
    k2 = k1 * 2
    k3 = (3 / m) * 3178.4
    k4 = k3 * 2
    k5 = 4 / (n - 1) * 22866.79
    k6 = 2 * k5
    k7 = 2 / (m - 1) * 23308.23
    k8 = 2 * k7
    # k_croix=(k6**2+k8**2)**(1/2)
    k_croix=10000 #je sais pas

    # longueurs au repos trouvees a partir du programme 5x3:
    l_repos_bord = 0.240-0.045233
    l_repos_coin = 0.240-0.100254
    l_repos_vertical = 4 * 1.052 / (n - 1)
    l_repos_horizontal = 2 * 1.0525 / (m - 1)
    l_croix=(l_repos_vertical**2+l_repos_horizontal**2)**0.5
    # #CALCUL DES RAIDEURS ET DES LONGUEURS AU REPOS

    # ressorts entre le cadre du trampoline et la toile : k1,k2,k3,k4
    k_bord = np.zeros(Nb_ressorts_cadre)
    l_bord_tab = np.zeros(Nb_ressorts_cadre)
    # cotes verticaux :
    k_bord[0:n], k_bord[n + m:2 * n + m] = k2, k2
    l_bord_tab[0:n], l_bord_tab[n + m:2 * n + m] = l_repos_bord, l_repos_bord
    # cotes horizontaux :
    k_bord[n:n + m], k_bord[2 * n + m:2 * n + 2 * m] = k4, k4
    l_bord_tab[n:n + m], l_bord_tab[2 * n + m:2 * n + 2 * m] = l_repos_bord, l_repos_bord
    # coins :
    k_bord[0], k_bord[n - 1], k_bord[n + m], k_bord[2 * n + m - 1] = k1, k1, k1, k1
    k_bord[n], k_bord[n + m - 1], k_bord[2 * n + m], k_bord[2 * (n + m) - 1] = k3, k3, k3, k3
    l_bord_tab[0], l_bord_tab[n - 1], l_bord_tab[n + m], l_bord_tab[2 * n + m - 1] = l_repos_coin, l_repos_coin, l_repos_coin, l_repos_coin
    l_bord_tab[n], l_bord_tab[n + m - 1], l_bord_tab[2 * n + m], l_bord_tab[2 * (n + m) - 1] = l_repos_coin, l_repos_coin, l_repos_coin, l_repos_coin

    # ressorts horizontaux internes a la toile : k5,k6
    k_horizontaux = k6 * np.ones(n * (m - 1))
    k_horizontaux[0:n * m - 1:n] = k5  # ressorts horizontaux du bord DE LA TOILE en bas
    k_horizontaux[n - 1:n * (m - 1):n] = k5  # ressorts horizontaux du bord DE LA TOILE en haut
    l_horizontal_tab =  l_repos_horizontal * np.ones(n * (m - 1))

    # ressorts verticaux internes a la toile : k7,k8
    k_verticaux = k8 * np.ones(m * (n - 1))
    k_verticaux[0:m * (n - 1):m] = k7  # ressorts verticaux du bord DE LA TOILE a droite
    k_verticaux[m - 1:n * m - 1:m] = k7  # ressorts verticaux du bord DE LA TOILE a gauche
    l_vertical_tab = l_repos_vertical * np.ones(m * (n - 1))

    #ressorts obliques internes a la toile :
    l_repos_croix=l_croix*np.ones(Nb_ressorts_croix)
    k_croix_tab=k_croix*np.ones(Nb_ressorts_croix)

    k = np.append(k_horizontaux, k_verticaux)
    k = np.append(k_bord, k)
    l_repos = np.append(l_horizontal_tab, l_vertical_tab)
    l_repos = np.append(l_bord_tab, l_repos)

    # CALCUL DES MASSES : (pas pris en compte la masse ajoutee par lathlete)
    mcoin = 1.803  # masse d'un point se trouvant sur un coin de la toile
    mpetit = 0.5 * 5.695 / (m - 2)  # masse d'un point se trouvant sur le petit cote de la toile
    mgrand = 0.5 * 9.707 / (n - 2)  # masse d'un point se trouvant sur le grand cote de la toile
    mmilieu = 3 * 0.650 / ((n - 2) * (m - 2))  # masse d'un point se trouvant au milieu de la toile

    M = mmilieu * np.ones(n * m)  # on initialise toutes les masses a celle du centre
    M[0], M[n - 1], M[n * (m - 1)], M[n * m - 1] = mcoin, mcoin, mcoin, mcoin
    M[n:n * (m - 1):n] = mpetit  # masses du cote bas
    M[2 * n - 1:n * m - 1:n] = mpetit  # masses du cote haut
    M[1:n - 1] = mgrand  # masse du cote droit
    M[n * (m - 1) + 1:n * m - 1] = mgrand  # masse du cote gauche
    M[int((m*n-1)/2)]+=Masse_centre

    ###################################################################################################################

    Pt = Optimisation_toile(Masse_centre,Pt,Pt_ancrage, k, M, l_repos, Pos_repos,k_croix_tab,l_repos_croix)

    Spring_bout_1, Spring_bout_2= Spring_bouts(Pt, Pt_ancrage)
    Spring_bout_1_repos, Spring_bout_2_repos = Spring_bouts(Pos_repos, Pt_ancrage)
    Spring_bout_croix_1, Spring_bout_croix_2 = Spring_bouts_croix(Pt)
    Spring_bout_croix_1_repos, Spring_bout_croix_2_repos = Spring_bouts_croix(Pos_repos)

    Spring_bout_1, Spring_bout_2 = Spring_bout_1.T,Spring_bout_2.T
    Spring_bout_1_repos,Spring_bout_2_repos=Spring_bout_1_repos.T,Spring_bout_2_repos.T
    Spring_bout_croix_1, Spring_bout_croix_2= Spring_bout_croix_1.T, Spring_bout_croix_2.T
    Spring_bout_croix_1_repos, Spring_bout_croix_2_repos = Spring_bout_croix_1_repos.T, Spring_bout_croix_2_repos.T

    Vect_unit_dir_F = (Spring_bout_2 - Spring_bout_1) / np.linalg.norm(Spring_bout_2 - Spring_bout_1)
    F_spring = np.zeros((3, Nb_ressorts))
    for ispring in range(Nb_ressorts):
        F_spring[:, ispring] = Vect_unit_dir_F[ispring, :] * k[ispring] * (
                np.linalg.norm(Spring_bout_2[ispring, :] - Spring_bout_1[ispring, :]) - l_repos[ispring])

    Vect_unit_dir_F_croix = (Spring_bout_croix_2 - Spring_bout_croix_1) /cas.norm_fro(Spring_bout_croix_2 - Spring_bout_croix_1)
    F_spring_croix = np.zeros((3, Nb_ressorts_croix))
    for ispring in range(Nb_ressorts_croix):
        F_spring_croix[:, ispring] = Vect_unit_dir_F_croix[ispring, :] * k_croix_tab[ispring] * (
                cas.norm_fro(Spring_bout_croix_2[ispring, :] - Spring_bout_croix_1[ispring, :]) - l_repos_croix[ispring])

    F_masses = np.zeros((3, n*m))
    F_masses[2, :] = - M * 9.81

    # return Pos_repos,Pt_ancrage, M, k, l_repos, k_croix_tab,l_repos_croix, F_spring, F_spring_croix,F_masses
    return M,F_spring, F_spring_croix,F_masses


def Calc_Spring_bouts (Pos_repos,Pt,Pt_ancrage) :
    Spring_bout_1, Spring_bout_2 = Spring_bouts(Pt, Pt_ancrage)
    Spring_bout_1_repos, Spring_bout_2_repos = Spring_bouts(Pos_repos, Pt_ancrage)
    Spring_bout_croix_1, Spring_bout_croix_2 = Spring_bouts_croix(Pt)
    Spring_bout_croix_1_repos, Spring_bout_croix_2_repos = Spring_bouts_croix(Pos_repos)

    Spring_bout_1, Spring_bout_2 = Spring_bout_1.T, Spring_bout_2.T
    Spring_bout_1_repos, Spring_bout_2_repos = Spring_bout_1_repos.T, Spring_bout_2_repos.T
    Spring_bout_croix_1, Spring_bout_croix_2 = Spring_bout_croix_1.T, Spring_bout_croix_2.T
    Spring_bout_croix_1_repos, Spring_bout_croix_2_repos = Spring_bout_croix_1_repos.T, Spring_bout_croix_2_repos.T

    return Spring_bout_1, Spring_bout_2, Spring_bout_1_repos, Spring_bout_2_repos, Spring_bout_croix_1, Spring_bout_croix_2, Spring_bout_croix_1_repos, Spring_bout_croix_2_repos


def Force_point(Pt,Masse_centre) : #--> resultante des forces en chaque point a un instant donne
    # Pos_repos, Pt_ancrage, M, k, l_repos, k_croix_tab,l_repos_croix, F_spring, F_spring_croix, F_masses = Force_calc(Pt,Masse_centre)

    M,F_spring, F_spring_croix, F_masses = Force_calc(Pt, Masse_centre)
    Pt_ancrage=Points_ancrage()

    # Spring_bout_1, Spring_bout_2, Spring_bout_1_repos, Spring_bout_2_repos, Spring_bout_croix_1, Spring_bout_croix_2, Spring_bout_croix_1_repos, Spring_bout_croix_2_repos = Calc_Spring_bouts (Pos_repos,Pt,Pt_ancrage)

    #forces de masse :
    F_masses = np.zeros((3, n * m))
    F_masses[2, :] = - M * 9.81

    #forces elastiques projetees sur x,y,z
    F_spring_points = np.zeros((3, n * m))

    # - points des coin de la toile : VERIFIE CEST OK
    F_spring_points[:,0]=F_spring[:,0]+\
                         F_spring[:,Nb_ressorts_cadre-1]+\
                         F_spring[:,Nb_ressorts_cadre]+ \
                         F_spring[:,Nb_ressorts_cadre+Nb_ressorts_horz] +\
                         F_spring_croix[:,0]# en bas a droite : premier ressort du cadre + dernier ressort du cadre + premiers ressorts horz, vert et croix
    F_spring_points[:, n-1] = F_spring[:,n-1] +\
                              F_spring[:, n] + \
                              F_spring[:, Nb_ressorts_cadre + n - 1] + \
                              F_spring[:, Nb_ressorts_cadre + Nb_ressorts_horz + Nb_ressorts_vert-m] + \
                              F_spring_croix[:,2*(n-1)-1]  # en haut a droite
    F_spring_points[:, (m-1)*n] = F_spring[:, 2*n+m-1] +\
                                  F_spring[:, 2*n+m] + \
                                  F_spring[:, Nb_ressorts_cadre + Nb_ressorts_horz-n] + \
                                  F_spring[:, Nb_ressorts_cadre + Nb_ressorts_horz + m-1] + \
                                  F_spring_croix[:, Nb_ressorts_croix - 2*(n-1) +1]  # en bas a gauche
    F_spring_points[:, m* n-1] = F_spring[:, n + m - 1] + \
                                 F_spring[:, n + m ] + \
                                 F_spring[:, Nb_ressorts_cadre + Nb_ressorts_horz-1] + \
                                 F_spring[:, Nb_ressorts-1] + \
                                 F_spring_croix[:, Nb_ressorts_croix-2]  # en haut a gauche

    # - points du bord de la toile> Pour lordre des termes de la somme, on part du ressort cadre puis sens trigo
            # - cote droit VERIFIE CEST OK
    for i in range (1,n-1):
        F_spring_points[:, i] = F_spring[:, i] + \
                                F_spring[:,Nb_ressorts_cadre + Nb_ressorts_horz + m * i] + \
                                F_spring_croix[:, 2 * (i - 1) + 1] + \
                                F_spring[:, Nb_ressorts_cadre + i] + \
                                F_spring_croix[:, 2 * (i - 1)+2] + \
                                F_spring[:,Nb_ressorts_cadre + Nb_ressorts_horz + m * (i - 1)]
            # - cote gauche VERIFIE CEST OK
    j=0
    for i in range((m-1)*n+1, m*n-1):
        F_spring_points[:,i]=F_spring[:,Nb_ressorts_cadre - m - (2+j)] + \
                             F_spring[:,Nb_ressorts_cadre+Nb_ressorts_horz+(j+1)*m-1]+ \
                             F_spring_croix[:,Nb_ressorts_croix-2*n+1+2*(j+2)]+\
                             F_spring[:,Nb_ressorts_cadre+Nb_ressorts_horz-n+j+1]+\
                             F_spring_croix[:,Nb_ressorts_croix-2*n+2*(j+1)]+\
                             F_spring[:,Nb_ressorts_cadre+Nb_ressorts_horz+(j+2)*m-1]
        j+=1

            # - cote haut VERIFIE CEST OK
    j=0
    for i in range (2*n-1,(m-1)*n,n) :
        F_spring_points[:, i]= F_spring[:, n+1+j] + \
                               F_spring[:, Nb_ressorts_cadre + i] + \
                               F_spring_croix[:,(j+2)*(n-1)*2-1]+\
                               F_spring[:,Nb_ressorts_cadre + Nb_ressorts_horz + (Nb_ressorts_vert+1) - (m-j)] +\
                               F_spring_croix[:,(j+1)*(n-1)*2-2]+\
                               F_spring[:, Nb_ressorts_cadre + i-n]
        j+=1
            # - cote bas VERIFIE CEST OK
    j=0
    for i in range (n,(m-2)*n+1,n) :
        F_spring_points[:, i] = F_spring[:, Nb_ressorts_cadre-(2+j)] + \
                                F_spring[:, Nb_ressorts_cadre + n*j]+\
                                F_spring_croix[:,1+2*(n-1)*j]+\
                                F_spring[:,Nb_ressorts_cadre+Nb_ressorts_horz+j+1]+\
                                F_spring[:,2*(n-1)*(j+1)]+\
                                F_spring[:, Nb_ressorts_cadre + n*(j+1)]
        j+=1

    #Points du centre de la toile (tous les points qui ne sont pas en contact avec le cadre)
    #on fait une colonne puis on passe a la colonne de gauche etc
    #dans lordre de la somme : ressort horizontal de droite puis sens trigo
    for j in range (1,m-1):
        for i in range (1,n-1) :
            F_spring_points[:,j*n+i]=F_spring[:,Nb_ressorts_cadre+(j-1)*n+i] + \
                                     F_spring_croix[:,2*j*(n-1) - 2*n + 3 + 2*i]+\
                                     F_spring[:,Nb_ressorts_cadre+Nb_ressorts_horz + m*i + j]+\
                                     F_spring_croix[:,j*2*(n-1) + i*2]+\
                                     F_spring[:, Nb_ressorts_cadre + j * n + i]+\
                                     F_spring_croix[:,j*2*(n-1) + i*2 -1]+\
                                     F_spring[:,Nb_ressorts_cadre+Nb_ressorts_horz + m*(i-1) + j]+\
                                     F_spring_croix[:,j*2*(n-1) -2*n + 2*i]
    F_point=F_masses+F_spring_points
    return F_point

###############################################################################
#INIT DYNAMIQUE :
dt = 1
T_total = 5

Pt_ancrage=Points_ancrage()
Pos_repos = np.zeros((3, n*m))
for i in range(n*m):
    Pos_repos[:, i] = np.array([l_toile - (i // n) * dl, -L_toile + (i % n) * dL, 0])

M,F_spring, F_spring_croix, F_masses = Force_calc(Pos_repos,Masse_centre)
F_point=Force_point(Pos_repos,Masse_centre)


Pt = np.zeros((int(T_total / dt), 3, np.shape(Pos_repos)[1]))
vitesse = np.zeros((int(T_total / dt), 3, np.shape(Pos_repos)[1]))
accel = np.zeros((int(T_total / dt), 3, np.shape(Pos_repos)[1]))


# F_point=Force_point(Pt,Masse_centre)

Pt[0, :, :] = Pos_repos[:, :]
print(Pt[0, :, :].shape)
Spring_bout_1, Spring_bout_2, Spring_bout_1_repos, Spring_bout_2_repos, Spring_bout_croix_1, Spring_bout_croix_2,Spring_bout_croix_1_repos, Spring_bout_croix_2_repos = Calc_Spring_bouts (Pos_repos,Pt[0,:,:],Pt_ancrage)

#######################################################################################################################
#AFFICHAGE ETAT INITIAL

fig = plt.figure()

ax = fig.add_subplot(1,int(T_total / dt),1, projection='3d')
ax.set_box_aspect([1.1, 1.8, 1])
ax.plot(Pos_repos[0, :], Pos_repos[1, :], Pos_repos[2, :], '.b')
ax.plot(Pt_ancrage[0, :], Pt_ancrage[1, :], Pt_ancrage[2, :], '.k')
# point du milieu  :
ax.plot(Pos_repos[0, int((n * m - 1) / 2)], Pos_repos[1, int((n * m - 1) / 2)], Pos_repos[2, int((n * m - 1) / 2)], '.y')

# ressorts entre le cadre et la toile :
# for j in range (2*(m+n)):
for j in range(Nb_ressorts):
    #pqs tres elegant mais cest le seul moyen pour que ca fonctionne
    a = []
    a = np.append(a, Spring_bout_1_repos[j, 0])
    a = np.append(a, Spring_bout_2_repos[j, 0])

    b = []
    b = np.append(b, Spring_bout_1_repos[j, 1])
    b = np.append(b, Spring_bout_2_repos[j, 1])

    c = []
    c = np.append(c, Spring_bout_1_repos[j, 2])
    c = np.append(c, Spring_bout_2_repos[j, 2])

    ax.plot3D(a, b, c, '-r', linewidth=1)

for j in range(Nb_ressorts_croix):
    # pqs tres elegant mais cest le seul moyen pour que ca fonctionne
    a = []
    a = np.append(a, Spring_bout_croix_1_repos[j, 0])
    a = np.append(a, Spring_bout_croix_2_repos[j, 0])

    b = []
    b = np.append(b, Spring_bout_croix_1_repos[j, 1])
    b = np.append(b, Spring_bout_croix_2_repos[j, 1])

    c = []
    c = np.append(c, Spring_bout_croix_1_repos[j, 2])
    c = np.append(c, Spring_bout_croix_2_repos[j, 2])

    ax.plot3D(a, b, c, '-g',linewidth=1)
plt.title('t=0 ')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')

####################################################################################################################

for t in range(1, int(T_total / dt)):
    F_point = Force_point(Pt,Masse_centre)
    accel[t, :, :] = F_point[:,:] / M[:]
    vitesse[t, :, :] = dt * accel[t, :, :] + vitesse[t - 1, :, :]
    Pt[t, :, :] = dt * vitesse[t, :, :] + Pt[t - 1, :, :]

    #######################################################################################################################
    #affichage a chaque instant
    ax = fig.add_subplot(1, int(T_total / dt), t+1, projection='3d')
    ax.set_box_aspect([1.1, 1.8, 1])
    ax.plot(Pt[t,0, :], Pt[t,1, :], Pt[t,2, :], '.b')
    ax.plot(Pt_ancrage[0, :], Pt_ancrage[1, :], Pt_ancrage[2, :], '.k')
    # point du milieu :
    ax.plot(Pt[t,0, int((n * m - 1) / 2)], Pt[t,1, int((n * m - 1) / 2)], Pt[t,2, int((n * m - 1) / 2)], '.y')

    ind_milieu_droite = int((m * n - 1) / 2 - n)
    ind_milieu_haut = int((m * n - 1) / 2 + 1)
    ind_milieu_gauche = int((m * n - 1) / 2 + n)
    ind_milieu_bas = int((m * n - 1) / 2 - 1)

    ind_milieu_bas_droite = int((m * n - 1) / 2 - n - 1)
    ind_milieu_haut_droite = int((m * n - 1) / 2 - n + 1)
    ind_milieu_haut_gauche = int((m * n - 1) / 2 + n + 1)
    ind_milieu_bas_gauche = int((m * n - 1) / 2 + n - 1)

    ax.plot(Pt[t,0, ind_milieu_droite], Pt[t,1, ind_milieu_droite], Pt[t,2, ind_milieu_droite], '.y')
    ax.plot(Pt[t,0, ind_milieu_haut], Pt[t,1, ind_milieu_haut], Pt[t,2, ind_milieu_haut], '.y')
    ax.plot(Pt[t,0, ind_milieu_gauche], Pt[t,1, ind_milieu_gauche], Pt[t,2, ind_milieu_gauche], '.y')
    ax.plot(Pt[t,0, ind_milieu_bas], Pt[t,1, ind_milieu_bas], Pt[t,2, ind_milieu_bas], '.y')

    ax.plot(Pt[t,0, ind_milieu_bas_droite], Pt[t,1, ind_milieu_bas_droite], Pt[t,2, ind_milieu_bas_droite], '.y')
    ax.plot(Pt[t,0, ind_milieu_haut_droite], Pt[t,1, ind_milieu_haut_droite], Pt[t,2, ind_milieu_haut_droite], '.y')
    ax.plot(Pt[t,0, ind_milieu_haut_gauche], Pt[t,1, ind_milieu_haut_gauche], Pt[t,2, ind_milieu_haut_gauche], '.y')
    ax.plot(Pt[t,0, ind_milieu_bas_gauche], Pt[t,1, ind_milieu_bas_gauche], Pt[t,2, ind_milieu_bas_gauche], '.y')

    for j in range(Nb_ressorts):
        # pqs tres elegant mais cest le seul moyen pour que ca fonctionne
        a = []
        a = np.append(a, Spring_bout_1[j, 0])
        a = np.append(a, Spring_bout_2[j, 0])

        b = []
        b = np.append(b, Spring_bout_1[j, 1])
        b = np.append(b, Spring_bout_2[j, 1])

        c = []
        c = np.append(c, Spring_bout_1[j, 2])
        c = np.append(c, Spring_bout_2[j, 2])

        ax.plot3D(a, b, c, '-r', linewidth=1)

    for j in range(Nb_ressorts_croix):
        # pqs tres elegant mais cest le seul moyen pour que ca fonctionne
        a = []
        a = np.append(a, Spring_bout_croix_1[j, 0])
        a = np.append(a, Spring_bout_croix_2[j, 0])


        b = []
        b = np.append(b, Spring_bout_croix_1[j, 1])
        b = np.append(b, Spring_bout_croix_2[j, 1])

        c = []
        c = np.append(c, Spring_bout_croix_1[j, 2])
        c = np.append(c, Spring_bout_croix_2[j, 2])

        ax.plot3D(a, b, c, '-g', linewidth=1)

    plt.title(
        'temps = ' +str(t))
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
plt.show()