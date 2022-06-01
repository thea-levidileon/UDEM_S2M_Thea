import numpy as np
# import numpy.matlib
# import sys
#import ipopt
# sys.path.append('/home/user/anaconda3/envs/qpoases/lib')
# sys.path.append('/home/user/anaconda3/envs/qpoases/include/casadi/build/lib')
# import casadi as cas
# from casadi import MX, SX, sqrt
# import biorbd
import matplotlib.pyplot as plt
# from os.path import dirname, join as pjoin
# import scipy.io as sio
# from IPython import embed
# from mpl_toolkits import mplot3d

n=15 #nombre de mailles sur le grand cote
m=9#nombre de mailles sur le petit cote
Masse_centre=150

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


def Points_ancrage_repos (Nb_increments) :
    # ancrage :
    Pt_ancrage = np.zeros((Nb_increments,Nb_ressorts_cadre,3))
    # cote droit :
    for i in range(n):
        Pt_ancrage[:,i,:] = np.array([l_cadre, -L_toile + i * dL, 0])
    # cote haut :
    for j in range(n, n + m):
        Pt_ancrage[:,j,:] = np.array([l_toile - (j - n) * dl, L_cadre, 0])
    # cote gauche :
    for k in range(n + m, 2 * n + m):
        Pt_ancrage[:,k,:] = np.array([-l_cadre, L_toile - (k - m - n) * dL, 0])
    # cote bas :
    for h in range(2 * n + m, Nb_ressorts_cadre):
        Pt_ancrage[:,h,:] = np.array([-l_toile + (h - 2 * n - m) * dl, -L_cadre, 0])

    #repos
    Pos_repos = np.zeros((Nb_increments,n*m, 3))
    for i in range(n*m):
        Pos_repos[:,i,:] = np.array([l_toile - (i // n) * dl, -L_toile + (i % n) * dL, 0])

    return Pt_ancrage,Pos_repos

def Param () :
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
    k_croix = 10000  # je sais pas

    # longueurs au repos trouvees a partir du programme 5x3:
    l_repos_bord = 0.240 - 0.045233
    l_repos_coin = 0.240 - 0.100254
    l_repos_vertical = 4 * 1.052 / (n - 1)
    l_repos_horizontal = 2 * 1.0525 / (m - 1)
    l_croix = (l_repos_vertical ** 2 + l_repos_horizontal ** 2) ** 0.5
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
    l_bord_tab[0], l_bord_tab[n - 1], l_bord_tab[n + m], l_bord_tab[
        2 * n + m - 1] = l_repos_coin, l_repos_coin, l_repos_coin, l_repos_coin
    l_bord_tab[n], l_bord_tab[n + m - 1], l_bord_tab[2 * n + m], l_bord_tab[
        2 * (n + m) - 1] = l_repos_coin, l_repos_coin, l_repos_coin, l_repos_coin

    # ressorts horizontaux internes a la toile : k5,k6
    k_horizontaux = k6 * np.ones(n * (m - 1))
    k_horizontaux[0:n * m - 1:n] = k5  # ressorts horizontaux du bord DE LA TOILE en bas
    k_horizontaux[n - 1:n * (m - 1):n] = k5  # ressorts horizontaux du bord DE LA TOILE en haut
    l_horizontal_tab = l_repos_horizontal * np.ones(n * (m - 1))

    # ressorts verticaux internes a la toile : k7,k8
    k_verticaux = k8 * np.ones(m * (n - 1))
    k_verticaux[0:m * (n - 1):m] = k7  # ressorts verticaux du bord DE LA TOILE a droite
    k_verticaux[m - 1:n * m - 1:m] = k7  # ressorts verticaux du bord DE LA TOILE a gauche
    l_vertical_tab = l_repos_vertical * np.ones(m * (n - 1))

    # ressorts obliques internes a la toile :
    l_repos_croix = l_croix * np.ones(Nb_ressorts_croix)
    k_croix_tab = k_croix * np.ones(Nb_ressorts_croix)

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
    M[int((m * n - 1) / 2)] += Masse_centre

    return k, l_repos, M, k_croix_tab, l_repos_croix

def Spring_bouts(Pt,Pt_ancrage,time,Nb_increments):
    # Definition des ressorts (position, taille)
    Spring_bout_1=np.zeros((Nb_increments,Nb_ressorts,3))

    # RESSORTS ENTRE LE CADRE ET LA TOILE
    for i in range(0, Nb_ressorts_cadre):
        Spring_bout_1[time,i,:] = Pt_ancrage[time-1,i, :]

    # RESSORTS HORIZONTAUX : il y en a n*(m-1)
    for i in range(Nb_ressorts_horz):
        Spring_bout_1[time,Nb_ressorts_cadre + i,:] = Pt[time-1,i,:]

    # RESSORTS VERTICAUX : il y en a m*(n-1)
    k=0
    for i in range(n - 1):
        for j in range(m):
            Spring_bout_1[time,Nb_ressorts_cadre+Nb_ressorts_horz+k, :] = Pt[time-1,i + n * j,:]
            k+=1
####################################################################################################################
    Spring_bout_2=np.zeros((Nb_increments,Nb_ressorts,3))

    # RESSORTS ENTRE LE CADRE ET LA TOILE
    for i in range(0, n): # points droite du bord de la toile
        Spring_bout_2[time,i,:] = Pt[time-1,i,:]

    k=0
    for i in range(n - 1, m * n, n): # points hauts du bord de la toile
        Spring_bout_2[time, n+k, :] = Pt[time-1, i, :]
        k+=1

    k=0
    for i in range(m*n-1,n * (m - 1)-1, -1): # points gauche du bord de la toile
        Spring_bout_2[time, n + m + k, :] = Pt[time-1, i, :]
        k+=1

    k=0
    for i in range(n * (m - 1), -1, -n): # points bas du bord de la toile
        Spring_bout_2[time, 2*n + m + k, :] = Pt[time-1, i, :]
        k+=1

    # RESSORTS HORIZONTAUX : il y en a n*(m-1)
    k=0
    for i in range(n, n * m):
        Spring_bout_2[time,Nb_ressorts_cadre + k,:] = Pt[time-1,i,:]
        k+=1

    # RESSORTS VERTICAUX : il y en a m*(n-1)
    k=0
    for i in range(1, n):
        for j in range(m):
            Spring_bout_2[time,Nb_ressorts_cadre + Nb_ressorts_horz + k,:] = Pt[time-1,i + n * j,:]
            k+=1

    return (Spring_bout_1,Spring_bout_2)

def Spring_bouts_croix(Pt,time,Nb_increments):
    #RESSORTS OBLIQUES : il n'y en a pas entre le cadre et la toile
    Spring_bout_croix_1=np.zeros((Nb_increments,Nb_ressorts_croix,3))

    #Pour spring_bout_1 on prend uniquement les points de droite des ressorts obliques
    k=0
    for i in range ((m-1)*n):
        Spring_bout_croix_1[time,k,:]=Pt[time-1,i,:]
        k += 1
        #a part le premier et le dernier de chaque colonne, chaque point est relie a deux ressorts obliques
        if (i+1)%n!=0 and i%n!=0 :
            Spring_bout_croix_1[time, k, :] = Pt[time-1, i, :]
            k+=1

    Spring_bout_croix_2=np.zeros((Nb_increments,Nb_ressorts_croix,3))
    #Pour spring_bout_2 on prend uniquement les points de gauche des ressorts obliques
    #pour chaue carre on commence par le point en haut a gauche, puis en bas a gauche
    #cetait un peu complique mais ca marche, faut pas le changer
    j=1
    k = 0
    while j<m:
        for i in range (j*n,(j+1)*n-2,2):
            Spring_bout_croix_2[time,k,:] = Pt[time-1,i + 1,:]
            Spring_bout_croix_2[time,k+1,:] = Pt[time-1,i,:]
            Spring_bout_croix_2[time,k+2,:] = Pt[time-1,i+ 2,:]
            Spring_bout_croix_2[time,k+3,:] = Pt[time-1,i + 1,:]
            k += 4
        j+=1

    return Spring_bout_croix_1,Spring_bout_croix_2


def ancien_Force_calc () :
    a=1

# def Force_calc(Pt,Masse_centre,time,Nb_increments):
#
#     Pt_ancrage,Pos_repos=Points_ancrage_repos(time,Nb_increments)
#     k, l_repos, M, k_croix_tab, l_repos_croix = Param()
#
#     Spring_bout_1, Spring_bout_2, Spring_bout_1_repos, Spring_bout_2_repos, Spring_bout_croix_1, Spring_bout_croix_2, \
#     Spring_bout_croix_1_repos, Spring_bout_croix_2_repos = Calc_Spring_bouts (Pos_repos,Pt,Pt_ancrage,time,Nb_increments)
#
#     F_spring = np.zeros((Nb_increments,3, Nb_ressorts))
#     Vect_unit_dir_F[time] = (Spring_bout_2[time] - Spring_bout_1[time]) / np.linalg.norm(Spring_bout_2[time] - Spring_bout_1[time])
#     for ispring in range(Nb_ressorts):
#         F_spring[time,:, ispring] = Vect_unit_dir_F[time,ispring, :] * k[ispring] * (
#                 np.linalg.norm(Spring_bout_2[time,ispring, :] - Spring_bout_1[time,ispring, :]) - l_repos[time,ispring])
#
#     F_spring_croix = np.zeros((Nb_increments,3, Nb_ressorts_croix))
#     Vect_unit_dir_F_croix = (Spring_bout_croix_2[time] - Spring_bout_croix_1[time]) /np.linalg.norm(Spring_bout_croix_2[time] - Spring_bout_croix_1[time])
#     for ispring in range(Nb_ressorts_croix):
#         F_spring_croix[time,:, ispring] = Vect_unit_dir_F_croix[time,ispring, :] * k_croix_tab[ispring] * (
#                 np.linalg.norm(Spring_bout_croix_2[time,ispring, :] - Spring_bout_croix_1[time,ispring, :]) - l_repos_croix[ispring])
#
#     F_masses = np.zeros((3, n*m))
#     F_masses[2, :] = - M * 9.81
#
#   return M,F_spring, F_spring_croix,F_masses
    return a


def Force_calc(Spring_bout_1,Spring_bout_2,Spring_bout_croix_1,Spring_bout_croix_2,Masse_centre,time,Nb_increments):

    # Pt_ancrage,Pos_repos=Points_ancrage_repos(time,Nb_increments)
    k, l_repos, M, k_croix_tab, l_repos_croix = Param()

    # Spring_bout_1, Spring_bout_2, Spring_bout_1_repos, Spring_bout_2_repos, Spring_bout_croix_1, Spring_bout_croix_2, \
    # Spring_bout_croix_1_repos, Spring_bout_croix_2_repos = Calc_Spring_bouts (Pos_repos,Pt,Pt_ancrage,time,Nb_increments)

    F_spring = np.zeros((Nb_ressorts,3))
    Vect_unit_dir_F = (Spring_bout_2 - Spring_bout_1) / np.linalg.norm(Spring_bout_2 - Spring_bout_1)
    for ispring in range(Nb_ressorts):
        F_spring[ispring,:] = Vect_unit_dir_F[ispring, :] * k[ispring] * (
                np.linalg.norm(Spring_bout_2[ispring, :] - Spring_bout_1[ispring, :]) - l_repos[ispring])

    F_spring_croix = np.zeros((Nb_ressorts_croix,3))
    Vect_unit_dir_F_croix = (Spring_bout_croix_2 - Spring_bout_croix_1) /np.linalg.norm(Spring_bout_croix_2 - Spring_bout_croix_1)
    for ispring in range(Nb_ressorts_croix):
        F_spring_croix[ispring,:] = Vect_unit_dir_F_croix[ispring, :] * k_croix_tab[ispring] * (
                np.linalg.norm(Spring_bout_croix_2[ispring, :] - Spring_bout_croix_1[ispring, :]) - l_repos_croix[ispring])

    F_masses = np.zeros((n*m,3))
    F_masses[:,2] = - M * 9.81

    return M,F_spring, F_spring_croix,F_masses


def Force_point(F_spring,F_spring_croix,F_masses,time,Nb_increments) : #--> resultante des forces en chaque point a un instant donne

    # M,F_spring, F_spring_croix, F_masses = Force_calc(Pt, Masse_centre,time,Nb_increments)

    #forces elastiques
    F_spring_points = np.zeros((n*m,3))

    # - points des coin de la toile : VERIFIE CEST OK
    F_spring_points[0,:]=F_spring[0,:]+\
                         F_spring[Nb_ressorts_cadre-1,:]+\
                         F_spring[Nb_ressorts_cadre,:]+ \
                         F_spring[Nb_ressorts_cadre+Nb_ressorts_horz,:] +\
                         F_spring_croix[0,:]# en bas a droite : premier ressort du cadre + dernier ressort du cadre + premiers ressorts horz, vert et croix
    F_spring_points[n-1,:] = F_spring[n-1,:] +\
                              F_spring[n,:] + \
                              F_spring[ Nb_ressorts_cadre + n - 1,:] + \
                              F_spring[ Nb_ressorts_cadre + Nb_ressorts_horz + Nb_ressorts_vert-m,:] + \
                              F_spring_croix[2*(n-1)-1,:]  # en haut a droite
    F_spring_points[ (m-1)*n,:] = F_spring[ 2*n+m-1,:] +\
                                  F_spring[ 2*n+m,:] + \
                                  F_spring[ Nb_ressorts_cadre + Nb_ressorts_horz-n,:] + \
                                  F_spring[ Nb_ressorts_cadre + Nb_ressorts_horz + m-1,:] + \
                                  F_spring_croix[ Nb_ressorts_croix - 2*(n-1) +1,:]  # en bas a gauche
    F_spring_points[ m* n-1,:] = F_spring[ n + m - 1,:] + \
                                 F_spring[ n + m,: ] + \
                                 F_spring[ Nb_ressorts_cadre + Nb_ressorts_horz-1,:] + \
                                 F_spring[ Nb_ressorts-1,:] + \
                                 F_spring_croix[ Nb_ressorts_croix-2,:]  # en haut a gauche

    # - points du bord de la toile> Pour lordre des termes de la somme, on part du ressort cadre puis sens trigo
            # - cote droit VERIFIE CEST OK
    for i in range (1,n-1):
        F_spring_points[ i,:] = F_spring[ i,:] + \
                                F_spring[Nb_ressorts_cadre + Nb_ressorts_horz + m * i,:] + \
                                F_spring_croix[ 2 * (i - 1) + 1,:] + \
                                F_spring[ Nb_ressorts_cadre + i,:] + \
                                F_spring_croix[ 2 * (i - 1)+2,:] + \
                                F_spring[Nb_ressorts_cadre + Nb_ressorts_horz + m * (i - 1),:]
            # - cote gauche VERIFIE CEST OK
    j=0
    for i in range((m-1)*n+1, m*n-1):
        F_spring_points[i,:]=F_spring[Nb_ressorts_cadre - m - (2+j),:] + \
                             F_spring[Nb_ressorts_cadre+Nb_ressorts_horz+(j+1)*m-1,:]+ \
                             F_spring_croix[Nb_ressorts_croix-2*n+1+2*(j+2),:]+\
                             F_spring[Nb_ressorts_cadre+Nb_ressorts_horz-n+j+1,:]+\
                             F_spring_croix[Nb_ressorts_croix-2*n+2*(j+1),:]+\
                             F_spring[Nb_ressorts_cadre+Nb_ressorts_horz+(j+2)*m-1,:]
        j+=1

            # - cote haut VERIFIE CEST OK
    j=0
    for i in range (2*n-1,(m-1)*n,n) :
        F_spring_points[ i,:]= F_spring[ n+1+j,:] + \
                               F_spring[ Nb_ressorts_cadre + i,:] + \
                               F_spring_croix[(j+2)*(n-1)*2-1,:]+\
                               F_spring[Nb_ressorts_cadre + Nb_ressorts_horz + (Nb_ressorts_vert+1) - (m-j),:] +\
                               F_spring_croix[(j+1)*(n-1)*2-2,:]+\
                               F_spring[ Nb_ressorts_cadre + i-n,:]
        j+=1
            # - cote bas VERIFIE CEST OK
    j=0
    for i in range (n,(m-2)*n+1,n) :
        F_spring_points[ i,:] = F_spring[ Nb_ressorts_cadre-(2+j),:] + \
                                F_spring[ Nb_ressorts_cadre + n*j,:]+\
                                F_spring_croix[1+2*(n-1)*j,:]+\
                                F_spring[Nb_ressorts_cadre+Nb_ressorts_horz+j+1,:]+\
                                F_spring[2*(n-1)*(j+1),:]+\
                                F_spring[ Nb_ressorts_cadre + n*(j+1),:]
        j+=1

    #Points du centre de la toile (tous les points qui ne sont pas en contact avec le cadre)
    #on fait une colonne puis on passe a la colonne de gauche etc
    #dans lordre de la somme : ressort horizontal de droite puis sens trigo
    for j in range (1,m-1):
        for i in range (1,n-1) :
            F_spring_points[j*n+i,:]=F_spring[Nb_ressorts_cadre+(j-1)*n+i,:] + \
                                     F_spring_croix[2*j*(n-1) - 2*n + 3 + 2*i,:]+\
                                     F_spring[Nb_ressorts_cadre+Nb_ressorts_horz + m*i + j,:]+\
                                     F_spring_croix[j*2*(n-1) + i*2,:]+\
                                     F_spring[ Nb_ressorts_cadre + j * n + i,:]+\
                                     F_spring_croix[j*2*(n-1) + i*2 -1,:]+\
                                     F_spring[Nb_ressorts_cadre+Nb_ressorts_horz + m*(i-1) + j,:]+\
                                     F_spring_croix[j*2*(n-1) -2*n + 2*i,:]
    F_point=F_masses+F_spring_points
    return F_point

def Etat_initial(Pt_ancrage,Pos_repos,Nb_increments,fig) :
    # Pt_ancrage, Pos_repos = Points_ancrage_repos(Nb_increments)
    Pt[0, :, :] = Pos_repos[0, :, :]

    Spring_bout_1, Spring_bout_2 = Spring_bouts(Pt, Pt_ancrage, 1, Nb_increments)
    Spring_bout_croix_1, Spring_bout_croix_2 = Spring_bouts_croix(Pt, 1, Nb_increments)
    Spb1, Spb2 = Spring_bout_1[1, :, :], Spring_bout_2[1, :, :]
    Spbc1, Spbc2 = Spring_bout_croix_1[1, :, :], Spring_bout_croix_2[1, :, :]

    ax = fig.add_subplot(int(T_total / (2 * dt)), int(T_total / (2 * dt)), 1, projection='3d')
    # ax = fig.add_subplot(1,1, 1, projection='3d')
    ax.set_box_aspect([1.1, 1.8, 1])
    ax.plot(Pt[0, :, 0], Pt[0, :, 1], Pt[0, :, 2], '.b')
    ax.plot(Pt_ancrage[0, :, 0], Pt_ancrage[0, :, 1], Pt_ancrage[0, :, 2], '.k')

    for j in range(Nb_ressorts):
        # pqs tres elegant mais cest le seul moyen pour que ca fonctionne
        a = []
        a = np.append(a, Spb1[j, 0])
        a = np.append(a, Spb2[j, 0])

        b = []
        b = np.append(b, Spb1[j, 1])
        b = np.append(b, Spb2[j, 1])

        c = []
        c = np.append(c, Spb1[j, 2])
        c = np.append(c, Spb2[j, 2])

        ax.plot3D(a, b, c, '-r', linewidth=1)

    for j in range(Nb_ressorts_croix):
        # pqs tres elegant mais cest le seul moyen pour que ca fonctionne
        a = []
        a = np.append(a, Spbc1[j, 0])
        a = np.append(a, Spbc2[j, 0])

        b = []
        b = np.append(b, Spbc1[j, 1])
        b = np.append(b, Spbc2[j, 1])

        c = []
        c = np.append(c, Spbc1[j, 2])
        c = np.append(c, Spbc2[j, 2])

        ax.plot3D(a, b, c, '-g', linewidth=1)

    plt.title('temps = ' + str(0))
    ax.axes.set_xlim3d(left=-2, right=2)
    ax.axes.set_ylim3d(bottom=-3, top=3)
    ax.axes.set_zlim3d(bottom=-1000, top=1000)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    #plt.show()
    return Spb1,Spb2,Spbc1,Spbc2

def Affichage(Pt,Pt_ancrage,Spb1,Spb2,Spbc1,Spbc2,time,fig) :
    ax = fig.add_subplot(int(T_total / (2 * dt)), int(T_total / (2 * dt)), time + 1, projection='3d')
    ax.set_box_aspect([1.1, 1.8, 1])
    ax.plot(Pt[time, :, 0], Pt[time, :, 1], Pt[time, :, 2], '.b')
    ax.plot(Pt_ancrage[time, :, 0], Pt_ancrage[time, :, 1], Pt_ancrage[time, :, 2], '.k')

    for j in range(Nb_ressorts):
        a = []
        a = np.append(a, Spb1[j, 0])
        a = np.append(a, Spb2[j, 0])

        b = []
        b = np.append(b, Spb1[j, 1])
        b = np.append(b, Spb2[j, 1])

        c = []
        c = np.append(c, Spb1[j, 2])
        c = np.append(c, Spb2[j, 2])

        ax.plot3D(a, b, c, '-r', linewidth=1)

    for j in range(Nb_ressorts_croix):
        a = []
        a = np.append(a, Spbc1[j, 0])
        a = np.append(a, Spbc2[j, 0])

        b = []
        b = np.append(b, Spbc1[j, 1])
        b = np.append(b, Spbc2[j, 1])

        c = []
        c = np.append(c, Spbc1[j, 2])
        c = np.append(c, Spbc2[j, 2])

        ax.plot3D(a, b, c, '-g', linewidth=1)

    plt.title('temps = ' + str(time) + ' s')
    # ax.axes.set_xlim3d(left=-2, right=2)
    # ax.axes.set_ylim3d(bottom=-3, top=3)
    # ax.axes.set_zlim3d(bottom=-1000, top=1000)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    # plt.show()

###############################################################################

#PARAMETRES POUR LA DYNAMIQUE :
dt = 1
T_total = 4
Nb_increments=int(T_total/dt)

Pt = np.zeros((Nb_increments, n*m,3))
vitesse = np.zeros((Nb_increments, n*m,3))
accel = np.zeros((Nb_increments, n*m,3))

Pt_ancrage, Pos_repos = Points_ancrage_repos(Nb_increments)
#######################################################################################################################

#BOUCLE TEMPORELLE
fig = plt.figure()

Spb1,Spb2,Spbc1,Spbc2 = Etat_initial(Pt_ancrage,Pos_repos,Nb_increments,fig) #--> actualise Pt[0,:,:] et fait l'affichage

for time in range(1, Nb_increments):
    M,F_spring, F_spring_croix, F_masses = Force_calc(Spb1,Spb2,Spbc1,Spbc2,
                                                    Masse_centre, time,Nb_increments)  # utilise Param()
    F_point = Force_point(F_spring, F_spring_croix, F_masses, time, Nb_increments)
    for index in range (n*m) :
        accel[time, index, :] = F_point[index,:] / M[index]
    vitesse[time, :, :] = dt * accel[time, :, :] + vitesse[time - 1, :, :]
    Pt[time, :, :] = dt * vitesse[time, :, :] + Pt[time - 1, :, :]

    Spring_bout_1, Spring_bout_2 = Spring_bouts(Pt, Pt_ancrage, time, Nb_increments)
    Spring_bout_croix_1, Spring_bout_croix_2 = Spring_bouts_croix(Pt, time, Nb_increments)
    Spb1, Spb2 = Spring_bout_1[time, :, :], Spring_bout_2[time, :, :]
    Spbc1, Spbc2 = Spring_bout_croix_1[time, :, :], Spring_bout_croix_2[time, :, :]

    Affichage(Pt, Pt_ancrage, Spb1, Spb2, Spbc1, Spbc2, time, fig)


plt.show()

