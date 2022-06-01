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


L=2.134
L_ressort=0.35
dL=2*L/(n-1)
l=1.07
l_ressort=0.38
dl=2*l/(m-1)


#repos :
Pos_repos = np.zeros((3, n*m + 8))
for i in range(n*m) :
	Pos_repos[:,i]=np.array([l-(i//n)*dl, -L+(i%n)*dL,0])
Pos_repos[:,n*m]=np.array([dl/2,0,0])
Pos_repos[:,n*m+1]=np.array([dl/2,dL/2,0])
Pos_repos[:,n*m+2]=np.array([0,dL/2,0])
Pos_repos[:,n*m+3]=np.array([-dl/2,dL/2,0])
Pos_repos[:,n*m+4]=np.array([-dl/2,0,0])
Pos_repos[:,n*m+5]=np.array([-dl/2,-dL/2,0])
Pos_repos[:,n*m+6]=np.array([0,-dL/2,0])
Pos_repos[:,n*m+7]=np.array([dl/2,-dL/2,0])


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

#points voisins du centre :
Pt_centre=np.zeros((3,5+8))
Pt_centre[:,0]=Pos_repos[:,int((m*n-1)/2)]
Pt_centre[:,1]=Pos_repos[:,int((m*n-1)/2-n)]
Pt_centre[:,2]=Pos_repos[:,int((m*n-1)/2+1)]
Pt_centre[:,3]=Pos_repos[:,int((m*n-1)/2+n)]
Pt_centre[:,4]=Pos_repos[:,int((m*n-1)/2-1)]
for j in range (5,13) :
    Pt_centre[:, j] = Pos_repos[:, n*m-1 + j-4]

#Affichage :


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

#trace des positions au repos :
plt.plot(Pos_repos[0,:],Pos_repos[1,:],'.',label='positions au repos')

#trace des points d'ancrage :
plt.plot(Pt_ancrage[0,:],Pt_ancrage[1,:],'.',color='red',label='points dancrage')

#trace des points autour du point central
plt.plot(Pt_centre[0,:],Pt_centre[1,:],'.',color='orange',label='points du centre')

plt.legend()
plt.title('Maillage '+str(m)+'x'+str(n) + ', avec '+str((m-2)*(n-2))+ '+ 8 capteurs a utiliser')
plt.axis('equal')
plt.show()