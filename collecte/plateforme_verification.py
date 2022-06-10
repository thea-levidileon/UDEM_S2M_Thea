
import numpy as np
from ezc3d import c3d
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import seaborn as sns

# trial_name = 'statique_vide'
# trial_name = 'statique_D8'
trial_name = 'statique_D11'
# trial_name = 'statiqueLeft_D8'
# trial_name = 'statiqueLeft_D11'
c3d_experimental = c3d(trial_name+'.c3d')

numero=2 #0:Fx, 1:Fy, 2:Fz, 3:Mx, 4:My, 5:Mz

#on garde seulement les pins qui ont un nom
force_labels=c3d_experimental['parameters']['ANALOG']['LABELS']['value']
indices_supp=[]
for i in range (len(force_labels)) :
    if 'Raw' in force_labels[i] :
        indices_supp=np.append(indices_supp,i)
ind_stop=int(indices_supp[0])

#affichage des donnees des plateformes
ana=c3d_experimental['data']['analogs'][0,0:ind_stop,:] #F1x,F1y,F1z,M1x,M1y,M1z,.....M4z size 24
for i in range (4) :
    index=numero+6*i
    longueur = np.size(ana[numero+6*i])
    print(str(c3d_experimental['parameters']['ANALOG']['LABELS']['value'][0:ind_stop][index]) + ' = ' + str(
        ana[index, int(longueur / 2)]))








    # parametres des plateformes
    # Nb_platform=c3d_experimental['parameters']['FORCE_PLATFORM']['USED']['value'][0]
    # platform_type=c3d_experimental['parameters']['FORCE_PLATFORM']['TYPE']['value']
    # platform_zero=c3d_experimental['parameters']['FORCE_PLATFORM']['ZERO']['value'] #pas normal ca doit etre des valeurs non nulles
    # platform_corners=c3d_experimental['parameters']['FORCE_PLATFORM']['CORNERS']['value'] #position des coins de chaeu plateforme : [xyz,quel coin,quelle plateforme] --> pb plateforme 1
    # platform_origin=c3d_experimental['parameters']['FORCE_PLATFORM']['ORIGIN']['value']
    # platform_channel=c3d_experimental['parameters']['FORCE_PLATFORM']['CHANNEL']['value'] #dit quel chanel contient quelle donnee
    # platform_calmatrix=c3d_experimental['parameters']['FORCE_PLATFORM']['CAL_MATRIX']['value'] #matrice de calibration : il n'y en a pas
