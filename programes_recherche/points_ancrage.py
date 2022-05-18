
import numpy as np
import numpy.matlib
import sys
# sys.path.append('/home/user/anaconda3/envs/qpoases/lib')
# sys.path.append('/home/user/anaconda3/envs/qpoases/include/casadi/build/lib')
import casadi as cas
# from casadi import MX, SX, sqrt
# import biorbd
import matplotlib.pyplot as plt
# from os.path import dirname, join as pjoin
# import scipy.io as sio
#from IPython import embed

n=10
m=4


L=2.134
L_ressort=0.35
dL=2*L/(n-1)
l=1.07
l_ressort=0.38
dl=2*l/(m-1)
    
Pt_ancrage = np.zeros((3, 2*(n+m)))
#cote droit :
for i in range (n):
    Pt_ancrage[:,i]=np.array([l+l_ressort,-L + i*dL,0])
#cote haut :
for j in range (n,n+m):
    Pt_ancrage[:,j]=np.array([l-(j-n)*dl,L+L_ressort,0])
    
#cote gauche :	
for k in range (n+m,2*n+m) :
    Pt_ancrage[:,k]=np.array([-l-l_ressort,L-(k-m-n)*dL,0])
    
#cote bas :   
for h in range (2*n+m,2*n+2*m) :
    Pt_ancrage[:,h]=np.array([-l+(h-2*n-m)*dl,-L-L_ressort,0])
    	


#trace du contour du cadre :
plt.plot([l+l_ressort,l+l_ressort],[-L-L_ressort,L+L_ressort],color='red')
plt.plot([-l-l_ressort,-l-l_ressort],[-L-L_ressort,L+L_ressort],color='red')
plt.plot([-l-l_ressort,l+l_ressort],[+L+L_ressort,L+L_ressort],color='red')
plt.plot([-l-l_ressort,l+l_ressort],[-L-L_ressort,-L-L_ressort],color='red')

#trace des points d'ancrage :
plt.plot(Pt_ancrage[0,:],Pt_ancrage[1,:],'o')



plt.show()

