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
# from matplotlib.patches import Rectangle
import matplotlib.patches as pat
import matplotlib.axes as axes

n=15
m=9

# L=2.134
L_ressort=0.35 #pas verifie
dLbord=0.267
dLmilieu=0.1524
dL=0.305
# L=(n-3)*dL/2 + dLbord
L=6*dL + dLbord

# l=1.07
l_ressort=0.38 #pas verifie
dl1=0.246
dl2=0.234
dlbord=0.1905
dlmilieu=0.117
l= 2*dl1 + dl2 + dlbord

#repos :
Pos_repos = np.zeros((3, n*m + 8))

#premiere colonne (celle de droite)
Pos_repos[:,0]=np.array([l,-L,0])
for i in range(0,n-2) :
    Pos_repos[:,i+1]=np.array([l,-L + dLbord + i*dL,0])
Pos_repos[:,n-1]=np.array([l,L,0])

#reste des colonnes
largeur=[l,l-dlbord,l-(dlbord+dl1),l-(dlbord+dl1+dl2),0,-(l-(dlbord+dl1+dl2)),-(l-(dlbord+dl1)),-(l-dlbord),-l]
for j in range (1,m):
    for i in range (0,n) :
        Pos_repos[1:2,j*n + i]=Pos_repos[1:2,i] #meme z et y que sur la premiere colonne
        Pos_repos[0, j * n + i] = largeur[j]

#Les 8 points du milieu :
largeur_milieu=[dlmilieu,dlmilieu,0,-dlmilieu,-dlmilieu,-dlmilieu,0,dlmilieu]
longueur_milieu=[0,dLmilieu,dLmilieu,dLmilieu,0,-dLmilieu,-dLmilieu,-dLmilieu]
for i in range (8) :
    Pos_repos[:, n * m + i] = np.array([largeur_milieu[i], longueur_milieu[i], 0])


#ancrage :
Pt_ancrage = np.zeros((3, 2*(n+m)))
# cote droit :
for i in range(n):
    Pt_ancrage[1:2, i] = Pos_repos[1:2,i]
    Pt_ancrage[0,i]=l + l_ressort
# cote haut :
for j in range(n, n + m):
    Pt_ancrage[:, j] = np.array([largeur[j-n], L + L_ressort, 0])
# cote gauche :
for k in range(n + m, 2 * n + m):
    Pt_ancrage[1:2, k] = - Pos_repos[1:2, k-n-m]
    Pt_ancrage[0,k]= -l - l_ressort
# cote bas :
for h in range(2 * n + m, 2 * n + 2 * m):
    Pt_ancrage[:,h]= np.array([-largeur[h-2*n-m], -L - L_ressort, 0])

#Affichage :
fig=plt.figure()
ax = fig.add_subplot(1, 1, 1)


#trace du contour du cadre :
plt.plot([l+l_ressort,l+l_ressort],[-L-L_ressort,L+L_ressort],color='red')
plt.plot([-l-l_ressort,-l-l_ressort],[-L-L_ressort,L+L_ressort],color='red')
plt.plot([-l-l_ressort,l+l_ressort],[+L+L_ressort,L+L_ressort],color='red')
plt.plot([-l-l_ressort,l+l_ressort],[-L-L_ressort,-L-L_ressort],color='red',label='cadre du trampoline')

#trace du contour de la toile :
plt.plot([l,l],[-L,L],color='green')
plt.plot([-l,-l],[-L,L],color='green')
plt.plot([-l,l],[L,L],color='green')
plt.plot([-l,l],[-L,-L],color='green',label='contour de la toile')

rect=plt.Rectangle((-l,-L), 2*l, 2*L,color='green',alpha=0.3)
ax.add_patch(rect)


#trace des positions au repos :
plt.plot(Pos_repos[0,:],Pos_repos[1,:],'.',label='positions au repos')

#trace des points d'ancrage :
plt.plot(Pt_ancrage[0,:],Pt_ancrage[1,:],'.',color='red',label='points dancrage')


plt.legend()
plt.title('Maillage '+str(m)+'x'+str(n) + ', avec '+str((m-2)*(n-2))+ '+ 8 marqueurs a utiliser')
plt.axis('equal')
plt.show()