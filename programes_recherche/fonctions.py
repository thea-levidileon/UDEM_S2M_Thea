import numpy as np
import numpy.matlib
import sys
# sys.path.append('/home/user/anaconda3/envs/qpoases/lib')
# sys.path.append('/home/user/anaconda3/envs/qpoases/include/casadi/build/lib')
import casadi as cas
#from casadi import MX, SX, sqrt
# import biorbd
import matplotlib.pyplot as plt
# from os.path import dirname, join as pjoin
# import scipy.io as sio
#from IPython import embed


def Energie_func(Pt, k, M, l_repos,Spring_bout_1,Spring_bout_2,Nb_ressorts):
    Energie = cas.MX.zeros(1)  # on initialise l'energie a 0 (valeur scalaire)
    for i in range(Nb_ressorts):
        Energie += 0.5 * k[i] * (cas.norm_fro(Spring_bout_2[:, i] - Spring_bout_1[:, i]) - l_repos[i]) ** 2  # E+=0.5*k*deltaL**2
    Energie += cas.sum2(M.T * 9.81 * Pt[2, :])  # a conserver E=mgz
    # return Energie
    func = cas.Function('Energie_func', [Pt], [Energie]).expand()
    return func

def Spring_bouts(Pt,Pt_ancrage,n,m):
    # Definition des ressorts (position, taille)
    Spring_bout_1 = []

    # RESSORTS ENTRE LE CADRE ET LA TOILE
    for i in range(0, 2 * m + 2 * n):
        Spring_bout_1 = cas.horzcat(Spring_bout_1, Pt_ancrage[:, i])

    # RESSORTS HORIZONTAUX : il y en a n*(m-1)
    for i in range(n * (m - 1)):
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
    for i in range(n * (m - 1), m * n):
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

def ForceEquilib_func(Pt,k, M, l_repos,Spring_bout_1,Spring_bout_2,Nb_ressorts):
    Vect_unit_dir_F = (Spring_bout_2 - Spring_bout_1) / cas.norm_fro(Spring_bout_2 - Spring_bout_1)
    F_spring = cas.MX.zeros(3, Nb_ressorts)
    for ispring in range(Nb_ressorts):
        F_spring[:, ispring] = Vect_unit_dir_F[:, ispring] * k[ispring] * (
                cas.norm_fro(Spring_bout_2[:, ispring] - Spring_bout_1[:, ispring]) - l_repos[ispring])
    # return F_spring
    func = cas.Function('ForceEquilib_func', [Pt], [F_spring])  # .expand()
    return func

def ForceEquilib_centre_func(Pt, k, M, l_repos, Masse_centre,Nb_ressorts):
    #k_centre = k[13, 18, 19, 24] #a changer
    #l_repos_centre = l_repos[13, 18, 19, 24] #a changer

    k_centre = k[
        2 * (n + m) + (m * n - 1) / 2 - n, 2 * (n + m) + n * (m - 1) + m * (n - 1) / 2 + (m - 1) / 2, 2 * (n + m) + (
                    m * n - 1) / 2, 2 * (n + m) + n * (m - 1) + m * ((n - 1) / 2 - 1) + (m - 1) / 2]
    l_repos_centre = l_repos[
        2 * (n + m) + (m * n - 1) / 2 - n, 2 * (n + m) + n * (m - 1) + m * (n - 1) / 2 + (m - 1) / 2, 2 * (n + m) + (
                    m * n - 1) / 2, 2 * (n + m) + n * (m - 1) + m * ((n - 1) / 2 - 1) + (m - 1) / 2]

    Pt_centre = np.zeros((3, 5))
    Pt_centre[:, 0] = Pt[:, int((m * n - 1) / 2)]
    Pt_centre[:, 1] = Pt[:, int((m * n - 1) / 2 - n)]
    Pt_centre[:, 2] = Pt[:, int((m * n - 1) / 2 + 1)]
    Pt_centre[:, 3] = Pt[:, int((m * n - 1) / 2 + n)]
    Pt_centre[:, 4] = Pt[:, int((m * n - 1) / 2 - 1)]

    Spring_bout_1_centre = cas.horzcat(Pt_centre[:, 1],Pt_centre[:, 2],Pt_centre[:, 3],Pt_centre[:, 4])
    Spring_bout_2_centre = cas.horzcat(Pt_centre[:, 0], Pt_centre[:, 0], Pt_centre[:, 0], Pt_centre[:, 0])

    Vect_unit_dir_F = (Spring_bout_1_centre - Spring_bout_2_centre) / cas.norm_fro(Spring_bout_1_centre - Spring_bout_2_centre)
    F_spring = cas.MX.zeros(3, Nb_ressorts)
    for ispring in range(4):
        F_spring[:, ispring] = Vect_unit_dir_F[:, ispring] * k_centre[ispring] * (
                    cas.norm_fro(Spring_bout_1_centre[:, ispring] - Spring_bout_2_centre[:, ispring]) - l_repos_centre[ispring])

    Force_Masse = cas.MX.zeros(3)
    Force_Masse[2] = - (Masse_centre + M[int((m * n - 1) / 2)]) * 9.81

    Force_tot = cas.sum2(F_spring) + Force_Masse

    func = cas.Function('ForceEquilib_func', [Pt], [Force_tot])  #.expand()

    return func

def L_const_func(Pt,Spring_bout_1,Spring_bout_2):

    Vect = cas.MX.zeros(Nb_ressorts)
    for ispring in range(Nb_ressorts):
        Vect[ispring] = cas.norm_fro(Spring_bout_2[:, ispring] - Spring_bout_1[:, ispring])

    func = cas.Function('L_const_func', [Pt], [Vect])  #.expand()

    return func

def Optimisation_toile(Masse_centre, Pt_ancrage, k, M, l_repos, Pos_repos):

    Pt = cas.MX.sym('Pt', 3, 15)

    Pt_ancrage = cas.MX(Pt_ancrage)
    k = cas.MX(k)
    M = cas.MX(M)
    l_repos = cas.MX(l_repos)

    Energie = Energie_func(Pt, k, M, l_repos, Spring_bout_1, Spring_bout_2, Nb_ressorts)
    ForceEquilib_centre = ForceEquilib_centre_func(Pt, k, M, l_repos, Masse_centre, Nb_ressorts)
    ForceEquilib = ForceEquilib_func(Pt, k, M, l_repos, Spring_bout_1, Spring_bout_2, Nb_ressorts)
    L_const = L_const_func(Pt,Spring_bout_1,Spring_bout_2)


    w = []  # vecteur de variables
    w0 = []  # conditions initiales
    lbw = []
    ubw = []
    g = []  # contraintes
    lbg = []
    ubg = []
    Pt_post = cas.MX.zeros(3, n*m)
    for i in range(n*m):
        Pt_var = cas.MX.sym(f'Pt_{i}', 3)
        Pt_post[:, i] = Pt_var
        w += [Pt_var]
        w0 += [Pos_repos[:, i]]
        lbw += [Pos_repos[0, i] - 3, Pos_repos[1, i] - 3, Pos_repos[2, i] - 3]  # [np.zeros(3) - 0.00001]  #
        ubw += [Pos_repos[0, i] + 3, Pos_repos[1, i] + 3, Pos_repos[2, i] + 3]  # [np.zeros(3)+ 0.00001]  #
            # if i == 7:
            #     g += [Pt_var[2]]
            #     lbg += [np.ones(1)*-3]
            #     ubg += [np.zeros(1)]

    g += [ForceEquilib_centre(Pt_post)]
    lbg += [-10 * np.ones(3)]  # np.zeros(3) - 1e-10]
    ubg += [10 * np.ones(3)]  # np.zeros(3) + 1e-10]

    F_masses = cas.MX.zeros(3, n*m)
    F_masses[2, :] = - M * 9.81

    F_spring = ForceEquilib(Pt_post)

        # g += [L_const(Pt_post) - l_repos]
        # lbg += [np.zeros(38)]
        # ubg += [np.ones(38)*10]

    Eq_Obj = cas.sum1((F_spring[:, 0] + F_spring[:, 6] + F_spring[:, 11] + F_spring[:, 5] + F_masses[:, 0]) ** 2 + \
                          (F_spring[:, 1] + F_spring[:, 7] + F_spring[:, 12] + F_spring[:, 6] + F_masses[:, 1]) ** 2 + \
                          (F_spring[:, 2] + F_spring[:, 8] + F_spring[:, 13] + F_spring[:, 7] + F_masses[:, 2]) ** 2 + \
                          (F_spring[:, 3] + F_spring[:, 9] + F_spring[:, 14] + F_spring[:, 8] + F_masses[:, 3]) ** 2 + \
                          (F_spring[:, 4] + F_spring[:, 10] + F_spring[:, 15] + F_spring[:, 9] + F_masses[:, 4]) ** 2 + \
                          (F_spring[:, 11] + F_spring[:, 17] + F_spring[:, 22] + F_spring[:, 16] + F_masses[:,
                                                                                                   5]) ** 2 + \
                          (F_spring[:, 12] + F_spring[:, 18] + F_spring[:, 23] + F_spring[:, 17] + F_masses[:,
                                                                                                   6]) ** 2 + \
                          (F_spring[:, 13] + F_spring[:, 19] + F_spring[:, 24] + F_spring[:, 18] + F_masses[:,
                                                                                                   7]) ** 2 + \
                          (F_spring[:, 14] + F_spring[:, 20] + F_spring[:, 25] + F_spring[:, 19] + F_masses[:,
                                                                                                   8]) ** 2 + \
                          (F_spring[:, 15] + F_spring[:, 21] + F_spring[:, 26] + F_spring[:, 20] + F_masses[:,
                                                                                                   9]) ** 2 + \
                          (F_spring[:, 22] + F_spring[:, 28] + F_spring[:, 33] + F_spring[:, 27] + F_masses[:,
                                                                                                   10]) ** 2 + \
                          (F_spring[:, 23] + F_spring[:, 29] + F_spring[:, 34] + F_spring[:, 28] + F_masses[:,
                                                                                                   11]) ** 2 + \
                          (F_spring[:, 24] + F_spring[:, 30] + F_spring[:, 35] + F_spring[:, 29] + F_masses[:,
                                                                                                   12]) ** 2 + \
                          (F_spring[:, 25] + F_spring[:, 31] + F_spring[:, 36] + F_spring[:, 30] + F_masses[:,
                                                                                                   13]) ** 2 + \
                          (F_spring[:, 26] + F_spring[:, 32] + F_spring[:, 37] + F_spring[:, 31] + F_masses[:,
                                                                                                   14]) ** 2)

    obj = Energie(Pt_post)  # + Eq_Obj

        # prob = {'f': obj, 'x': cas.vertcat(*w), 'g': cas.vertcat(*g)}

        # opts = {'ipopt.max_iter':1000, 'ipopt.bound_push':1e-10, 'ipopt.bound_frac':1e-10}

        # solver = cas.nlpsol('solver', 'ipopt', prob) # opts 'ipopt'
        # sol = solver(x0=cas.vertcat(*w0), lbg=cas.vertcat(*lbg), ubg=cas.vertcat(*ubg), lbx=cas.vertcat(*lbw), ubx=cas.vertcat(*ubw),)

    qp = {'x': cas.vertcat(*w), 'f': obj, 'g': cas.vertcat(*g)}
    solver = cas.qpsol('solver', 'qpoases', qp)  # chercersur intenret doc casadi, chercher exeples ipopt

    sol = solver(x0=cas.vertcat(*w0), lbg=cas.vertcat(*lbg), ubg=cas.vertcat(*ubg), lbx=cas.vertcat(*lbw),
                     ubx=cas.vertcat(*ubw))

    Solution = sol['x']

    Pt = np.zeros((3, n*m))
    for j in range(n*m):
        Pt[:, j] = np.reshape(Solution[3 * j: 3 * j + 3], (3))

        # Output_IPOPT = solver.stats()
        # if Output_IPOPT['success'] == True:
        #     np.save(f'Position_massPoints/Pos_{Pos_y}y_{Pos_z}z', Pt)
        # else:
        #     np.save(f'Position_massPoints/Pos_{Pos_y}y_{Pos_z}z_PAS_CONVERGE', Pt)

    return Pt

#######################################################################################################################
#Parametres geometriques :
n=15
m=9
Nb_ressorts=2*n*m+n+m #nombre de ressorts total dans le modele


########################################################################################################################
#Parametres physiques :
#k trouves a partir du programme 5x3:
k1=(5/(n))*3266.68
k2=k1*2
k3=(3/(m))*3178.4
k4=k3*2
k5=4/(n-1)*22866.79
k6=2*k5
k7=2/(m-1)*23308.23
k8=2*k7

#longueurs au repos trouvees a partir du programme 5x3:
l_bord=0.240
l_coin=0.240
l_vertical=4*1.052/(n-1)
l_horizontal=2*1.0525/(m-1)

# #CALCUL DES RAIDEURS ET DES LONGUEURS AU REPOS

#ressorts entre le cadre du trampoline et la toile : k1,k2,k3,k4
k_bord=np.zeros(2*(m+n))
l_bord_tab=np.zeros(2*(m+n))
#cotes verticaux :
k_bord[0:n],k_bord[n+m:2*n+m]=k2,k2
l_bord_tab[0:n],l_bord_tab[n+m:2*n+m]=l_bord,l_bord
#cotes horizontaux :
k_bord[n:n+m],k_bord[2*n+m:2*n+2*m] = k4,k4
l_bord_tab[n:n+m],l_bord_tab[2*n+m:2*n+2*m] = l_bord,l_bord
#coins :
k_bord[0],k_bord[n-1],k_bord[n+m],k_bord[2*n+m-1]=k1,k1,k1,k1
k_bord[n],k_bord[n+m-1],k_bord[2*n+m],k_bord[2*(n+m)-1]=k3,k3,k3,k3
l_bord_tab[0],l_bord_tab[n-1],l_bord_tab[n+m],l_bord_tab[2*n+m-1]=l_coin,l_coin,l_coin,l_coin
l_bord_tab[n],l_bord_tab[n+m-1],l_bord_tab[2*n+m],l_bord_tab[2*(n+m)-1]=l_coin,l_coin,l_coin,l_coin

#ressorts horizontaux internes a la toile : k5,k6
k_horizontaux=k6*np.ones(n*(m-1))
# for i in range (n*(m-1)):
#     if i%n==0 or i%n==n-1 :
#         k_horizontaux[i]=k5
k_horizontaux[0:n*m-1:n]=k5 #ressorts horizontaux du bord DE LA TOILE en bas
k_horizontaux[n-1:n*(m-1):n]=k5 #ressorts horizontaux du bord DE LA TOILE en haut
l_horizontal_tab=l_horizontal*np.ones(n*(m-1))

#ressorts verticaux internes a la toile : k7,k8
k_verticaux=k8*np.ones(m*(n-1))
# for i in range (m*(n-1)):
#     if i%m==0 or i%m==m-1 :
#         k_verticaux[i]=k7
k_verticaux[0:m*(n-1):m]=k7 #ressorts verticaux du bord DE LA TOILE a droite
k_verticaux[m-1:n*m-1:m]=k7 #ressorts verticaux du bord DE LA TOILE a gauche
l_vertical_tab=l_vertical*np.ones(m*(n-1))

k_milieu=np.append(k_horizontaux,k_verticaux)
k=np.append(k_bord,k_milieu)

l_milieu=np.append(l_horizontal_tab,l_vertical_tab)
l_repos=np.append(l_bord_tab,l_milieu)

#CALCUL DES MASSES : (pas pris en compte la masse ajoutee par lathlete)
mcoin=1.803 #masse d'un point se trouvant sur un coin de la toile
mpetit= 0.5*5.695/(m-2)#masse d'un point se trouvant sur le petit cote de la toile
mgrand=0.5*9.707/(n-2) #masse d'un point se trouvant sur le grand cote de la toile
mcentre=3*0.650/((n-2)*(m-2)) #masse d'un point se trouvant au milieu de la toile

M=mcentre*np.ones(n*m) #on initialise toutes les masses a celle du centre
M[0],M[n-1],M[n*(m-1)],M[n*m-1]=mcoin,mcoin,mcoin,mcoin
M[n:n*(m-1):n]=mpetit #masses du cote bas
M[2*n-1:n*m-1:n]=mpetit #masses du cote haut
M[1:n-1]=mgrand #masse du cote droit
M[n*(m-1)+1:n*m-1]=mgrand #masse du cote gauche

######################################################################################################################
#definition des points :

L=2.134
L_ressort=0.35
dL=2*L/(n-1)
l=1.07
l_ressort=0.38
dl=2*l/(m-1)


#repos :
Pos_repos = np.zeros((3, n*m))
for i in range(n*m) :
	Pos_repos[:,i]=np.array([l-(i//n)*dl, -L+(i%n)*dL,0])


#ancrage :
Pt_ancrage = np.zeros((3, 2*(n+m)))
# cote droit :
for i in range(n):
    Pt_ancrage[:, i] = np.array([l + l_ressort, -L + i * dL, 0])
# cote haut :
for j in range(n, n + m):
    Pt_ancrage[:, j] = np.array([l - (j - n) * dl, L + L_ressort, 0])
# cote gauche :
for g in range(n + m, 2 * n + m):
    Pt_ancrage[:, g] = np.array([-l - l_ressort, L - (g - m - n) * dL, 0])
# cote bas :
for h in range(2 * n + m, 2 * n + 2 * m):
    Pt_ancrage[:, h] = np.array([-l + (h - 2 * n - m) * dl, -L - L_ressort, 0])

######################################################################################################################

###########################################################################################################################
#Resultat :

# Spring_bout_1,Spring_bout_2=Spring_bouts(Pos_repos,Pt_ancrage,n,m)
# energie = Energie(Pos_repos, k, M, l_repos,Spring_bout_1,Spring_bout_2,Nb_ressorts)
# print(energie)
# F_spring=ForceEquilib_func(Pos_repos, k, M, l_repos,Spring_bout_1,Spring_bout_2,Nb_ressorts)
# print(F_spring)

#Affichage :
