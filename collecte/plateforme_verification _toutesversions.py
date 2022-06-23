
import numpy as np
from ezc3d import c3d
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import seaborn as sns
import scipy as sc
from scipy import signal

action = 'acceleration'#'calibrage' #'dynamique' #'soustraction' #'comparaison' #'acceleration'
calibration = 1 #0 #1
calcul_rapport=1
numero_disque = 6
participant=1


Nb_disque=['vide','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11']
Poids_disque=[0,197.181,394.362,592.524, 789.705,937.836,1135.017,1333.179,1531.341,1728.522,1926.684, 2174.877]

vide_name='labeled_statiqueVide03'
statique_name='labeled_statique_' + Nb_disque[numero_disque]
trial_name= 'SautHaut_03'
# trial_name='p1_statique02'
# trial_name='p1_piedscolles_01'
# trial_name='labeled_p1_sauthaut_01'
calibrage_name = 'Test_Plateformes_4empilees05'

c3d_vide = c3d(vide_name+'.c3d')
c3d_statique = c3d(statique_name+'.c3d')
c3d_experimental = c3d(trial_name+'.c3d')
c3d_calibrage = c3d(calibrage_name+'.c3d')


def Param_platform(c3d_experimental) :
    # parametres des plateformes
    Nb_platform=c3d_experimental['parameters']['FORCE_PLATFORM']['USED']['value'][0]
    platform_type=c3d_experimental['parameters']['FORCE_PLATFORM']['TYPE']['value']
    platform_zero=c3d_experimental['parameters']['FORCE_PLATFORM']['ZERO']['value'] #pas normal ca doit etre des valeurs non nulles
    platform_corners=c3d_experimental['parameters']['FORCE_PLATFORM']['CORNERS']['value'] #position des coins de chaeu plateforme : [xyz,quel coin,quelle plateforme] --> pb plateforme 1
    platform_origin=c3d_experimental['parameters']['FORCE_PLATFORM']['ORIGIN']['value']
    platform_channel=c3d_experimental['parameters']['FORCE_PLATFORM']['CHANNEL']['value'] #dit quel chanel contient quelle donnee
    platform_calmatrix=c3d_experimental['parameters']['FORCE_PLATFORM']['CAL_MATRIX']['value'] #matrice de calibration : il n'y en a pas
    return (Nb_platform,platform_type,platform_zero,platform_corners,platform_origin,platform_channel,platform_calmatrix)

def Named_pins(c3d_experimental) :
    #on garde seulement les Raw pins
    force_labels=c3d_experimental['parameters']['ANALOG']['LABELS']['value']
    ind=[]
    for i in range (len(force_labels)) :
        if 'Raw' in force_labels[i] :
            ind_stop=i
            ind=np.append(ind,ind_stop)
    ind_stop=int(ind[0])
    return ind_stop

def matrices() :

    #M1,M2,M4 sont les matrices obtenues apres la calibration sur la plateforme 3
    M4_old = [[4.1390, 0.4018, -0.0456, -0.0280, -0.1606, -0.9432],
          [-0.7248, 2.7855, 0.1181, 0.0009, 0.1933, 0.9311],
          [0.0111, -0.0063, 0.8371, -0.004, 0.0014, 0.0120],
          [-4.4537, 3.1811, -0.2881, 0.9530, -1.2086, -4.5486],
          [-0.1021, 0.0049, -0.0019, -0.0139, 0.2170, -0.0016],
          [-0.1899, -0.0652, 0.0083, 0.0113, 0.0341, 0.1495]]  # c'est la bonne matrice pour PF4 (vieille S2M)

    M4=[[2.8658,0.0071,0.0547,0.0002,-0.0684,0.3753],
        [0.2035,3.0273,0.0106,0.0223,0.0232,-0.0245],
        [-0.0169,-0.005,0.8406,0.0003,-0.0007,0.02],
        [-0.2262,0.1468,-0.0247,0.0594,-0.0911,-0.3481],
        [0.009,0.017,-0.0042,-0.0167,0.2080,-0.034],
        [0.0073,-0.0025,0.0074,-0.0014,0.0008,0.4829]]

    M1_old = [[1.5289740, -0.0055372, 0.0000385, -0.0073467, 0.0008197, -0.7711510],
          [0.0145983, 1.5235620, 0.0000019, -0.0216603, -0.0003118, 0.7447910],
          [0.0041178, 0.0013416, 0.3924790, -0.0051392, -0.0000419, -0.0002830],
          [0.0077071, 0.0930658, 0.0254123, 1.4258270, 0.0021123, 0.8327660],
          [0.0059914, 0.0018655, -0.0010224, 0.0029237, 0.0917925, 0.0435810],
          [0.0017551, 0.0050230, -0.0020866, 0.0109248, -0.0002057, 0.1652510]]  # New_S2M, matrice 7273M pin 10 a 15

    M1=[[1.2579,0.0943,0.0065,-0.0474,-0.0527,0.0104],
        [0.2055,1.3526,-0.0099,0.037,0.0375,0.0411],
        [0.025,0.0184,0.3898,-0.0032,-0.0035,0.0601],
        [-0.0401,-0.0946,0.0101,0.0691,0.0709,-0.2894],
        [0.0737,0.0573,-0.0092,0.0201,0.0234,0.2787],
        [0.0819,0.0778,-0.0079,-0.0219,-0.0228,0.1781]]

    M2_old = [[1.5378351, 0.0014677, -0.0003717, 0.0005122, -0.0002392, -0.0012910],
          [-0.0046930, 1.5345052, 0.0033831, 0.0002683, -0.0001478, -0.0006520],
          [-0.0016725, -0.0021274, 0.3923144, 0.0001368, 0.0000125, -0.0001880],
          [0.0019166, -0.0004030, 0.0012135, 0.0898106, -0.0000327, 0.0017920],
          [0.0008414, -0.0004323, -0.0003952, -0.0005261, 0.0897077, 0.0004990],
          [-0.0013925, -0.0007727, -0.0015802, -0.0000313, -0.0003715,
           0.1960360]]  # INS-7382, matrice 7382M-1-1 pin 19 a 24

    M2=[[1.635,-0.0277,0.0167,-0.0052,0.0172,-0.0613],
        [-0.0305,1.5842,-0.0121,0.0013,0.0029,0.0693],
        [-0.0047,-0.0012,0.3923,0.0007,-0.0008,0.0007],
        [-0.1530,0.1225,-0.0034,0.0176,-0.0589,-0.1458],
        [-0.0548,0.0247,-0.0033,-0.0176,0.0752,-0.0362],
        [-0.0039,-0.0009,0.0088,-0.0006,0.0011,0.2114]]

    M3 = [[1.5438977, -0.0073794, 0.0022438, 0.0000509, -0.0002489, -0.0003690],
          [0.0104010, 1.5370655, 0.0013308, -0.0001873, 0.0000820, 0.0011140],
          [-0.0022109, -0.0003121, 0.3901421, 0.0001339, -0.0000556, -0.0002320],
          [0.0010562, -0.0010583, 0.0019478, 0.0900776, 0.0000017, 0.0014160],
          [0.0009278, -0.0021898, -0.0002942, 0.0002198, 0.0900668, 0.0003340],
          [-0.0003099, 0.0000449, 0.0035554, -0.0005566, -0.0002628,
           0.1962430]]  # INS-7383, matrice 7383M-1-5 pin 28 a 33

    # zeros donnes par Nexus
    zeros1 = np.array([1.0751899, 2.4828501, -0.1168980, 6.8177500, -3.0313399, -0.9456340])
    zeros2 = np.array([0., -2., -2., 0., 0., 0.])
    zeros3 = np.array([0.0307411, -5., -4., -0.0093422, -0.0079338, 0.0058189])
    zeros4 = np.array([-0.1032560, -3., -3., 0.2141770, 0.5169040, -0.3714130])

    # return M1_old,M2_old,M3,M4_old,zeros1,zeros2,zeros3,zeros4
    return M1, M2, M3, M4, zeros1, zeros2, zeros3, zeros4

def plateforme_calcul (c3d) :
    ind_stop = Named_pins(c3d)
    ana = c3d['data']['analogs'][0, ind_stop:, :]
    M1, M2, M3, M4, zeros1, zeros2, zeros3, zeros4 = matrices()

    platform1 = ana[0:6, :]
    platform2 = ana[6:12, :]
    platform3 = ana[12:18, :]
    platform4 = ana[18:24, :]

    if calibration == 1 :
    #calibration des donnees brutes :
        platform1 = np.matmul(M1, platform1) * 100 * 10
        platform2 = np.matmul(M2, platform2) * 200 * 10
        platform3 = np.matmul(M3, platform3) * 100 * 10
        platform4 = np.matmul(M4, platform4) * 25 * 10

    platform=np.array([platform1, platform2, platform3, platform4])

    return platform

def soustraction (c3d_statique,c3d_vide) :
    platform_statique=plateforme_calcul (c3d_statique)
    platform_vide = plateforme_calcul(c3d_vide)

    for j in range (6) :
        for i in range (4) :
            platform_statique[i,j,:]=platform_statique[i,j,:] - np.mean(platform_vide[i,j])

    return platform_statique

def comparaison (c3d_statique,c3d_vide) :
    platform_statique = plateforme_calcul(c3d_statique)
    platform_vide = plateforme_calcul(c3d_vide)
    return platform_statique, platform_vide

def dynamique(c3d_experimental) :
    platform_experimental = plateforme_calcul(c3d_experimental)
    longueur = np.size(platform_experimental[0,0])
    zero_variable=np.zeros((4,6))
    for i in range (6) :
        for j in range(4) :
            zero_variable[j,i]=np.mean(platform_experimental[j,i,0:100])
            platform_experimental[j, i,:] = platform_experimental[j, i,:] - zero_variable[j,i]*np.ones(longueur)
    return platform_experimental

def calibrage(c3d_calibrage) :
    ind = [i for i in range(18, 36)] + [j for j in range(45, 51)]
    ana=c3d_calibrage['data']['analogs'][0, ind, :]

    M1, M2, M3, M4, zeros1, zeros2, zeros3, zeros4 = matrices()
    platform1 = ana[0:6, :]
    platform2 = ana[6:12, :]
    platform3 = ana[12:18, :]
    platform4 = ana[18:24, :]

    if calibration == 1: #calibration des donnees brutes :
        platform1 = np.matmul(M1, platform1) * 100 * 10
        platform2 = np.matmul(M2, platform2) * 200 * 10
        platform3 = np.matmul(M3, platform3) * 100 * 10
        platform4 = np.matmul(M4, platform4) * 25 * 10

    platform_calibrage = np.array([platform1, platform2, platform3, platform4])
    longueur = np.size(platform_calibrage[0, 0])
    zero_variable = np.zeros((4, 6))
    for i in range(6):
        for j in range(4):
            zero_variable[j, i] = np.mean(platform_calibrage[j, i, 0:100])
            platform_calibrage[j, i, :] = platform_calibrage[j, i, :] - zero_variable[j, i] * np.ones(longueur)

    return platform_calibrage

def unlabeled(c3d_experimental) :
    #boucle pour supprimer les points non labeled
    labels=c3d_experimental['parameters']['POINT']['LABELS']['value']
    indices_supp=[]
    for i in range (len(labels)) :
        if '*' in labels[i] :
            indices_supp=np.append(indices_supp,i)
    ind_stop=int(indices_supp[0])
    #labels et points avec les points non labelles supprimes
    labels=c3d_experimental['parameters']['POINT']['LABELS']['value'][0:ind_stop]
    points=c3d_experimental['data']['points'][:3,:ind_stop,:]
    return ind_stop,labels,points

def acceleration_point_bas(ind_stop,labels,points) : #trouver le point qui descend le plus bas :

    #on cherche le z min de chaque marqueur (on en profite pour supprimer les nan)
    minimum_marqueur=[np.nanmin(points[2,i]) for i in range (ind_stop)]
    argmin_temps=[np.where((points[2,i])==minimum_marqueur[i]) for i in range (ind_stop)]

    #on cherche le z min total
    argmin_marqueur=np.argmin(minimum_marqueur)
    argmin=[argmin_marqueur,int(argmin_temps[argmin_marqueur][0])]
    label_min=labels[argmin[0]]

    print('Altitude minimale obtenue au marqueur ' + str(label_min) + ' au temps t=' + str(argmin[1]))

    point_bas = points[:, argmin[0]]*0.001 #on veut les coordonnees en metres !!
    T = np.size(points[0, 0])
    dt=0.002
    vitesse_point_bas=np.array([[(point_bas[j, i + 1] - point_bas[j, i])/dt for i in range(0, T - 1)] for j in range(3)])
    accel_point_bas=np.array([[(vitesse_point_bas[j, i + 1] - vitesse_point_bas[j, i])/dt for i in range(0, T - 2)] for j in range(3)])
    # accel_point_bas = np.array([[(point_bas[j, i + 2] - 2 * point_bas[j, i + 1] + point_bas[j, i])/dt**2 for i in range(0, T - 2)] for j in range(3)])

    return label_min, point_bas, vitesse_point_bas, accel_point_bas, T


fig = plt.figure()
fig.add_subplot(313)

if action == 'soustraction' :
    platform = soustraction (c3d_statique,c3d_vide) #platefroem statique a laquelle on a soustrait la valeur de la plateforme a vide
    longueur = np.size(platform[0,0])
    x = np.linspace(0, longueur, longueur)

    plt.subplot(311)
    plt.plot(x, platform[0,0, :], '-k', label='Fx PF1')
    plt.plot(x, platform[1,0, :], '-b', label='Fx PF2')
    plt.plot(x, platform[2,0, :], '-r', label='Fx PF3')
    plt.plot(x, platform[3,0, :], '-m', label='Fx PF4')
    plt.legend()
    plt.title('Essai statique ' + str(statique_name) + ', valeurs a vide soustraites')

    plt.subplot(312)
    plt.plot(x, platform[0,1, :], '-k', label='Fy  PF1')
    plt.plot(x, platform[1,1, :], '-b', label='Fy PF2')
    plt.plot(x, platform[2,1, :], '-r', label='Fy PF3')
    plt.plot(x, platform[3,1, :], '-m', label='Fy PF4')
    plt.legend()

    plt.subplot(313)
    plt.plot(x, platform[0,2, :], '-k', label='Fz PF1')
    plt.plot(x, platform[1,2, :], '-b', label='Fz PF2')
    plt.plot(x, platform[2,2, :], '-r', label='Fz PF3')
    plt.plot(x, platform[3,2, :], '-m', label='Fz PF4')
    plt.plot(x,(Poids_disque[numero_disque] / 4)*np.ones(longueur), '-g',label='Poids theorique')
    plt.legend()

if action == 'comparaison' :
    platform_statique,platform_vide = comparaison(c3d_statique,c3d_vide)
    longueur = np.minimum(np.size(platform_statique[0, 0]),np.size(platform_vide[0, 0]))
    x = np.linspace(0, longueur, longueur)

    plt.subplot(311)
    plt.plot(x, platform_statique[0,0, :longueur], '-k', label='Fx PF1')
    plt.plot(x, platform_statique[1,0, :longueur], '-b', label='Fx PF2')
    plt.plot(x, platform_statique[2,0, :longueur], '-r', label='Fx PF3')
    plt.plot(x, platform_statique[3,0, :longueur], '-m', label='Fx PF4')
    plt.plot(x, platform_vide[0, 0, :longueur], '.k', label='Fx PF1')
    plt.plot(x, platform_vide[1, 0, :longueur], '.b', label='Fx PF2')
    plt.plot(x, platform_vide[2, 0, :longueur], '.r', label='Fx PF3')
    plt.plot(x, platform_vide[3, 0, :longueur], '.m', label='Fx PF4')
    plt.legend()
    plt.title('Essai statique ' + str(statique_name) + ' compare avec l\'essai a vide ' + str(vide_name))

    plt.subplot(312)
    plt.plot(x, platform_statique[0,1, :longueur], '-k', label='Fy PF1')
    plt.plot(x, platform_statique[1,1, :longueur], '-b', label='Fy PF2')
    plt.plot(x, platform_statique[2,1, :longueur], '-r', label='Fy PF3')
    plt.plot(x, platform_statique[3,1, :longueur], '-m', label='Fy PF4')
    plt.plot(x, platform_vide[0, 1, :longueur], '.k', label='Fy PF1')
    plt.plot(x, platform_vide[1, 1, :longueur], '.b', label='Fy PF2')
    plt.plot(x, platform_vide[2, 1, :longueur], '.r', label='Fy PF3')
    plt.plot(x, platform_vide[3, 1, :longueur], '.m', label='Fy PF4')
    plt.legend()

    plt.subplot(313)
    plt.plot(x, platform_statique[0,2, :longueur], '-k', label='Fz PF1')
    plt.plot(x, platform_statique[1,2, :longueur], '-b', label='Fz PF2')
    plt.plot(x, platform_statique[2,2, :longueur], '-r', label='Fz PF3')
    plt.plot(x, platform_statique[3,2, :longueur], '-m', label='Fz PF4')
    plt.plot(x, platform_vide[0, 2, :longueur], '.k', label='Fz PF2')
    plt.plot(x, platform_vide[1, 2, :longueur], '.b', label='Fz PF2')
    plt.plot(x, platform_vide[2, 2, :longueur], '.r', label='Fz PF3')
    plt.plot(x, platform_vide[3, 2, :longueur], '.m', label='Fz PF4')
    if calibration==1 :
        plt.plot(x, (Poids_disque[numero_disque]/4)*np.ones(longueur), '-g', label='Poids theorique divise par 4 pour ' + str(numero_disque) + ' disques')
    plt.legend()

if action == 'dynamique' :
    platform = dynamique(c3d_experimental)
    longueur = np.size(platform[0, 0])
    x = np.linspace(0, longueur, longueur)
    if calcul_rapport==0 :
        plt.subplot(311)
        plt.plot(x, platform[0, 0, :], '-k', label='Fx calibre PF1')
        plt.plot(x, platform[1, 0, :], '-b', label='Fx calibre PF2')
        plt.plot(x, platform[2, 0, :], '-r', label='Fx calibre PF3')
        plt.plot(x, platform[3, 0, :], '-m', label='Fx calibre PF4')
        plt.legend()
        plt.title('Essai dynamique ' + str(trial_name))

        plt.subplot(312)
        plt.plot(x, platform[0, 1, :], '-k', label='Fy calibre PF1')
        plt.plot(x, platform[1, 1, :], '-b', label='Fy calibre PF2')
        plt.plot(x, platform[2, 1, :], '-r', label='Fy calibre PF3')
        plt.plot(x, platform[3, 1, :], '-m', label='Fy calibre PF4')
        plt.legend()

        rapport=2.92
        poids = 64.5 * 9.81
        # poids=87.2*9.81
        plt.subplot(313)
        plt.plot(x, platform[0, 2, :]*rapport, '-k', label='Fz PF1 multiplie par rapport de poids')
        plt.plot(x, platform[1, 2, :]*rapport, '-b', label='Fz PF2 multiplie par rapport de poids')
        plt.plot(x, platform[2, 2, :]*rapport, '-r', label='Fz PF3 multiplie par rapport de poids')
        plt.plot(x, platform[3, 2, :]*rapport, '-m', label='Fz PF4 multiplie par rapport de poids')
        plt.plot(x, (platform[0, 2, :] + platform[1, 2, :] + platform[2, 2, :] + platform[3, 2, :]) * rapport, '-g',
                 label='moyenne multpliee par le rapport de poids')
        plt.plot(x, poids * np.ones(longueur), '-y', label='poids theorique')
        plt.legend()

    else :
        poids = 87.2 * 9.81
        # poids = 64.5 * 9.81
        debut_plateau =10000
        fin_plateau=35000
        moyenne_statique = np.mean(platform[0, 2, debut_plateau:fin_plateau] + platform[1, 2, debut_plateau:fin_plateau] + platform[2, 2, debut_plateau:fin_plateau] + platform[3, 2, debut_plateau:fin_plateau])
        moyenne_plateforme = np.mean(platform[0, 2, 16000:28000])
        # ind = []
        # for i in range(longueur):
        #     if platform[0, 2, i] + platform[1, 2, i] + platform[2, 2, i] + platform[3, 2, i] > 10:
        #         ind = np.append(ind, i)
        # debut=int(ind[0])
        # fin=int(ind[-1])
        # moyenne_tot=np.mean(platform[0, 2, debut:fin] + platform[1, 2, debut:fin] + platform[2, 2, debut:fin] + platform[3, 2,debut:fin])

        # rapport=moyenne_tot/moyenne
        # rapport = poids / moyenne_statique
        rapport=2.92

        plt.subplot(311)
        plt.plot(x, platform[0, 2, :], '-k', label='Fz PF1 sans modif')
        plt.plot(x, platform[1, 2, :], '-b', label='Fz PF2 sans modif')
        plt.plot(x, platform[2, 2, :], '-r', label='Fz PF3 sans modif')
        plt.plot(x, platform[3, 2, :], '-m', label='Fz PF4 sans modif')
        plt.plot(x, platform[0, 2, :] + platform[1, 2, :] + platform[2, 2, :] + platform[3, 2, :], '-g',
                 label='somme sans modif')
        plt.plot(x, poids * np.ones(longueur), '-y', label='poids theorique')
        plt.legend()
        plt.title('Calcul du rapport pour arriver au bon poids : ' +str(rapport) + '. Erreur = ' + str(moyenne_statique*rapport - poids))

        plt.subplot(312)
        plt.plot(x, platform[0, 2, :], '-k', label='Fz PF1 sans modif')
        plt.plot(x, platform[1, 2, :], '-b', label='Fz PF2 sans modif')
        plt.plot(x, platform[2, 2, :], '-r', label='Fz PF3 sans modif')
        plt.plot(x, platform[3, 2, :], '-m', label='Fz PF4 sans modif')
        # plt.plot(x, moyenne_plateforme * np.ones(longueur), '-m',label='moyenne d\'une plateforme')
        plt.plot(x,(platform[0, 2, :]+platform[1, 2, :]+platform[2, 2, :]+platform[3, 2, :])*rapport,'-g',label='somme multpliee par le rapport de poids')
        plt.plot(x,poids*np.ones(longueur),'-y',label='poids theorique')
        plt.legend()

        plt.subplot(313)
        plt.plot(x, platform[0, 2, :]*rapport, '-k', label='Fz PF1 multiplie par rapport')
        plt.plot(x, platform[1, 2, :]*rapport, '-b', label='Fz PF2 multiplie par rapport')
        plt.plot(x, platform[2, 2, :]*rapport, '-r', label='Fz PF3 multiplie par rapport')
        plt.plot(x, platform[3, 2, :]*rapport, '-m', label='Fz PF4 multiplie par rapport')
        # plt.plot(x, moyenne_plateforme * np.ones(longueur), '-m')
        plt.plot(x, (platform[0, 2, :] + platform[1, 2, :] + platform[2, 2, :] + platform[3, 2, :])*rapport, '-g',
                 label='somme multpliee par le rapport de poids')
        plt.plot(x, moyenne_statique * np.ones(longueur)*rapport, '-g', label='moyenne de la somme x rapport')
        plt.plot(x, poids * np.ones(longueur), '-y', label='poids theorique')
        plt.legend()

if action == 'calibrage' :
    platform=calibrage(c3d_calibrage)
    longueur = np.size(platform[0, 0])
    x = np.linspace(0, longueur, longueur)

    plt.subplot(311)
    plt.plot(x, platform[0, 0, :], '-k', label='Fx calibre PF1')
    plt.plot(x, platform[1, 0, :], '-b', label='Fx calibre PF2')
    plt.plot(x, platform[2, 0, :], '-r', label='Fx calibre PF3')
    plt.plot(x, platform[3, 0, :], '-m', label='Fx calibre PF4')
    plt.legend()
    plt.title('Fichier de calibration apres calibration : ' + str(calibrage_name))

    plt.subplot(312)
    plt.plot(x, platform[0, 1, :], '-k', label='Fy calibre  PF1')
    plt.plot(x, platform[1, 1, :], '-b', label='Fy calibre PF2')
    plt.plot(x, platform[2, 1, :], '-r', label='Fy calibre PF3')
    plt.plot(x, platform[3, 1, :], '-m', label='Fy calibre PF4')
    plt.legend()

    plt.subplot(313)
    plt.plot(x, platform[0, 2, :], '-k', label='Fz calibre PF1')
    plt.plot(x, platform[1, 2, :], '-b', label='Fz calibre PF2')
    plt.plot(x, platform[2, 2, :], '-r', label='Fz calibre PF3')
    plt.plot(x, platform[3, 2, :], '-m', label='Fz calibre PF4')
    plt.legend()

if action == 'acceleration' :
    rapport = 2.92

    if participant == 1 :
        masse = 64.5
    if participant == 2:
        masse = 87.2

    ind_stop, labels, points = unlabeled(c3d_experimental)
    label_min, point_bas, vitesse_point_bas, accel_point_bas, T = acceleration_point_bas(ind_stop, labels, points)

    platform = dynamique(c3d_experimental)

    x = np.linspace(0, T - 2, T - 2)

    plt.subplot(211)
    plt.plot(x, masse * accel_point_bas[2], '-c')
    plt.title('Acceleration verticale du point ' + str(label_min) + ' au cours de l\'essai ' + str(trial_name))

    plt.subplot(212)
    #on filtre l'accceleration :
    sos = signal.butter(10, [1,10], 'bandpass', fs=500, output='sos')
    accel_point_bas = signal.sosfilt(sos, accel_point_bas)
    plt.plot(x, masse*accel_point_bas[2], '-y', label='masse fois acceleration filtre')
    plt.plot(x,
             [(platform[0, 2, i * 4] + platform[1, 2, i * 4] + platform[2, 2, i * 4] + platform[3, 2, i * 4]) * rapport
              for i in range(T - 2)], '-g', label='somme des Fz des PF multipliee par le rapport de poids')
    plt.legend()

plt.show()



