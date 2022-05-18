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

n=15
m=9

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


#Affichage :
print('Masse totale : ' + str(np.sum(M)))
print('Longueur toile : '+ str(np.sum(l_vertical_tab)/m))
print('Largeur toile : '+ str(np.sum(l_horizontal_tab)/n))

plt.subplot(3,1,1)
xk=[i for i in range(np.size(k_bord) +np.size(k_horizontaux)+np.size(k_verticaux))]
plt.plot(xk,k,'-x')
plt.ylabel('Raideur du ressort (Nm)')
plt.xlabel('Indice du ressort')
plt.title('Raideurs des ressorts pour la maille ' + str(n) + 'x'+str(m)+' : ressorts dancrage puis ressorts horizontaux de la toile puis ressorts verticaux de la toile')

plt.subplot(3,1,2)
xl=[i for i in range(np.size(l_bord_tab) +np.size(l_horizontal_tab)+np.size(l_vertical_tab))]
plt.plot(xl,l_repos,'-x')
plt.ylabel('Longueur au repos du ressort (m)')
plt.xlabel('Indice du ressort')
plt.title('Longueur au repos des ressorts pour la maille ' + str(n) + 'x'+str(m)+' : ressorts dancrage puis ressorts horizontaux de la toile puis ressorts verticaux de la toile')

plt.subplot(3,1,3)
xm=[i for i in range(np.size(M))]
plt.plot(xm,M,'-x')
plt.ylabel('Masse (kg)')
plt.xlabel('Indice du point de masse')
plt.title('Masses localisees pour la maille ' + str(n) + 'x'+str(m)+' : lordre est le meme que pour les poisitions de repos' )

plt.show()