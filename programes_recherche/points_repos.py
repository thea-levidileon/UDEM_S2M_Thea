
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
from IPython import embed

n=5
m=3


L=2.134
L_ressort=0.35
dL=2*L/(n-1)
l=1.07
l_ressort=0.38
dl=2*l/(m-1)

Pos_repos = np.zeros((3, n*m))
for i in range(n*m) :
	Pos_repos[:,i]=np.array([l-(i//n)*dl, -L+(i%n)*dL,0])

#trace du contour de la toile :
plt.plot([l,l],[-L,L],color='green')
plt.plot([-l,-l],[-L,L],color='green')
plt.plot([-l,l],[L,L],color='green')
plt.plot([-l,l],[-L,-L],color='green',label='contour de la toile')

#trace des positions au repos :
plt.plot(Pos_repos[0,:],Pos_repos[1,:],'.',label='positions au repos')

plt.show()
