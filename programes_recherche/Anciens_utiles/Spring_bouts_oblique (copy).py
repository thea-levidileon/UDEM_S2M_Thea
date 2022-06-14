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
Nb_ressorts=2*n*m+n+m #nombre de ressorts total dans le modele
Nb_ressorts_croix=2*(m-1)*(n-1)


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
for k in range(n + m, 2 * n + m):
    Pt_ancrage[:, k] = np.array([-l - l_ressort, L - (k - m - n) * dL, 0])
# cote bas :
for h in range(2 * n + m, 2 * n + 2 * m):
    Pt_ancrage[:, h] = np.array([-l + (h - 2 * n - m) * dl, -L - L_ressort, 0])

###########################################################

Spring_bout_1=[]

#RESSORTS ENTRE LE CADRE ET LA TOILE
#tous les points dancrage
for i in range(0,2*m+2*n):
    Spring_bout_1=cas.horzcat(Spring_bout_1,Pt_ancrage[:,i])

#RESSORTS HORIZONTAUX : il y en a n*(m-1)
#pour spring_bout_1, on prend seulement les points de droite de chaque ressort horizontal
for i in range (n*(m-1)) :
    Spring_bout_1 = cas.horzcat(Spring_bout_1, Pos_repos[:, i])

#RESSORTS VERTICAUX : il y en a m*(n-1)
#pour spring_bout_1, on prend seulement les points du bas de chaque ressort vertical (de droite a gauche puis on remonte)
for i in range (n-1):
    for j in range (m) :
        Spring_bout_1 = cas.horzcat(Spring_bout_1, Pos_repos[:, i+n*j])

#####################################################

Spring_bout_2=[]

#RESSORTS ENTRE LE CADRE ET LA TOILE
#tous les points situes au bord de la toile
for i in range (0,n):
    Spring_bout_2=cas.horzcat(Spring_bout_2,Pos_repos[:,i]) #points droite du bord de la toile
for i in range (n-1,m*n,n):
    Spring_bout_2=cas.horzcat(Spring_bout_2,Pos_repos[:,i]) #points hauts du bord de la toile
for i in range (m*n-1,n * (m - 1)-1, -1):
    Spring_bout_2=cas.horzcat(Spring_bout_2,Pos_repos[:,i]) #points gauche du bord de la toile
for i in range(n*(m-1),-1,-n):
    Spring_bout_2 = cas.horzcat(Spring_bout_2, Pos_repos[:,i]) #points bas du bord de la toile

#RESSORTS HORIZONTAUX : il y en a n*(m-1)
#pour spring_bout_2, on prend seulment les points de gauche de chaque ressort horizontal
for i in range (n,n*m) :
    Spring_bout_2 = cas.horzcat(Spring_bout_2, Pos_repos[:, i])

#RESSORTS VERTICAUX : il y en a m*(n-1)
#pour spring_bout_1, on prend seulement les points du haut de chaque ressort vertical (de droite a gauche puis on remonte)
for i in range (1,n):
    for j in range (m) :
        Spring_bout_2 = cas.horzcat(Spring_bout_2, Pos_repos[:, i+n*j])






#RESSORTS OBLIQUES : il n'y en a pas entre le cadre et la toile

Spring_bout_croix_1=[]
#Pour spring_bout_1 on prend uniquement les points de droite des ressorts obliques
for i in range ((m-1)*n):
    Spring_bout_croix_1=cas.horzcat(Spring_bout_croix_1,Pos_repos[:,i])
    #a part le premier et le dernier de chaque colonne, chaque point est relie a deux ressorts obliques
    if (i+1)%n!=0 and i%n!=0 :
        Spring_bout_croix_1 = cas.horzcat(Spring_bout_croix_1, Pos_repos[:, i])

Spring_bout_croix_2=[]
#Pour spring_bout_2 on prend uniquement les points de gauche des ressorts obliques
#pour chaue carre on commence par le point en haut a gauche, puis en bas a gauche
#cetait un peu complique mais ca marche, faut pas le changer
j=1
while j<m:
    for i in range (j*n,(j+1)*n-2,2):
        Spring_bout_croix_2 = cas.horzcat(Spring_bout_croix_2, Pos_repos[:, i + 1])
        Spring_bout_croix_2 = cas.horzcat(Spring_bout_croix_2, Pos_repos[:, i])
        Spring_bout_croix_2 = cas.horzcat(Spring_bout_croix_2, Pos_repos[:, i+ 2])
        Spring_bout_croix_2 = cas.horzcat(Spring_bout_croix_2, Pos_repos[:, i + 1])
    j+=1


#Affichage :

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([1.1, 1.8, 1])
ax.plot(Pos_repos[0, :], Pos_repos[1, :], Pos_repos[2, :], '.b')
ax.plot(Pt_ancrage[0, :], Pt_ancrage[1, :], Pt_ancrage[2, :], '.k')
# point du milieu :
ax.plot(Pos_repos[0, int((n * m - 1) / 2)], Pos_repos[1, int((n * m - 1) / 2)], Pos_repos[2, int((n * m - 1) / 2)], '.y')

for j in range(Nb_ressorts):
    #pas tres elegant mais cest le seul moyen pour que ca fonctionne
    a = []
    a = np.append(a, Spring_bout_1.T[j, 0])
    a = np.append(a, Spring_bout_2.T[j, 0])

    b = []
    b = np.append(b, Spring_bout_1.T[j, 1])
    b = np.append(b, Spring_bout_2.T[j, 1])

    c = []
    c = np.append(c, Spring_bout_1.T[j, 2])
    c = np.append(c, Spring_bout_2.T[j, 2])

    ax.plot3D(a, b, c, '-r',linewidth=1)

for j in range(Nb_ressorts_croix):
    #pas tres elegant mais cest le seul moyen pour que ca fonctionne
    a = []
    a = np.append(a, Spring_bout_croix_1.T[j, 0])
    a = np.append(a, Spring_bout_croix_2.T[j, 0])

    b = []
    b = np.append(b, Spring_bout_croix_1.T[j, 1])
    b = np.append(b, Spring_bout_croix_2.T[j, 1])

    c = []
    c = np.append(c, Spring_bout_croix_1.T[j, 2])
    c = np.append(c, Spring_bout_croix_2.T[j, 2])

    ax.plot3D(a, b, c, '-g',linewidth=1)



plt.title('Trampoline sans sollicitation exterieure, maille ' +str(m)+'x'+str(n))
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')
plt.show()

# #trace du contour du cadre :
# plt.plot([l+l_ressort,l+l_ressort],[-L-L_ressort,L+L_ressort],color='red')
# plt.plot([-l-l_ressort,-l-l_ressort],[-L-L_ressort,L+L_ressort],color='red')
# plt.plot([-l-l_ressort,l+l_ressort],[+L+L_ressort,L+L_ressort],color='red')
# plt.plot([-l-l_ressort,l+l_ressort],[-L-L_ressort,-L-L_ressort],color='red',label='cadre du trampoline')
#
# #trace du contour de la toile :
# plt.plot([l,l],[-L,L],color='green')
# plt.plot([-l,-l],[-L,L],color='green')
# plt.plot([-l,l],[L,L],color='green')
# plt.plot([-l,l],[-L,-L],color='green',label='contour de la toile')

# #trace des positions au repos :
# plt.plot(Pos_repos[0,:],Pos_repos[1,:],Pos_repos[2,:],'.',label='positions au repos')
#
# #trace des points d'ancrage :
# plt.plot(Pt_ancrage[0,:],Pt_ancrage[1,:],Pt_ancrage[2,:],'.',color='red',label='points dancrage')

# plt.legend()
# plt.title('Maillage '+str(m)+'x'+str(n) + ', avec '+str((m-2)*(n-2))+ ' capteurs a utiliser')
# plt.axis('equal')
# plt.show()