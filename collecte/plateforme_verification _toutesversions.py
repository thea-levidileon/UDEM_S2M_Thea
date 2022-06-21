
import numpy as np
from ezc3d import c3d
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import seaborn as sns

action = 'comparaison'#'calibrage' #'dynamique' #'soustraction' #'comparaison'
calibration = 1 #0 #1
numero_disque = 11


Nb_disque=['vide','D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11']
Poids_disque=[0,197.181,394.362,592.524, 789.705,937.836,1135.017,1333.179,1531.341,1728.522,1926.684, 2174.877]

vide_name='labeled_statiqueVide03'
statique_name='labeled_statique_left_' + Nb_disque[numero_disque]
trial_name= 'SautHaut_03'
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
    M4 = [[4.1390, 0.4018, -0.0456, -0.0280, -0.1606, -0.9432],
          [-0.7248, 2.7855, 0.1181, 0.0009, 0.1933, 0.9311],
          [0.0111, -0.0063, 0.8371, -0.004, 0.0014, 0.0120],
          [-4.4537, 3.1811, -0.2881, 0.9530, -1.2086, -4.5486],
          [-0.1021, 0.0049, -0.0019, -0.0139, 0.2170, -0.0016],
          [-0.1899, -0.0652, 0.0083, 0.0113, 0.0341, 0.1495]]  # c'est la bonne matrice pour PF4 (vieille S2M)

    M1 = [[1.5289740, -0.0055372, 0.0000385, -0.0073467, 0.0008197, -0.7711510],
          [0.0145983, 1.5235620, 0.0000019, -0.0216603, -0.0003118, 0.7447910],
          [0.0041178, 0.0013416, 0.3924790, -0.0051392, -0.0000419, -0.0002830],
          [0.0077071, 0.0930658, 0.0254123, 1.4258270, 0.0021123, 0.8327660],
          [0.0059914, 0.0018655, -0.0010224, 0.0029237, 0.0917925, 0.0435810],
          [0.0017551, 0.0050230, -0.0020866, 0.0109248, -0.0002057, 0.1652510]]  # New_S2M, matrice 7273M pin 10 a 15

    M2 = [[1.5378351, 0.0014677, -0.0003717, 0.0005122, -0.0002392, -0.0012910],
          [-0.0046930, 1.5345052, 0.0033831, 0.0002683, -0.0001478, -0.0006520],
          [-0.0016725, -0.0021274, 0.3923144, 0.0001368, 0.0000125, -0.0001880],
          [0.0019166, -0.0004030, 0.0012135, 0.0898106, -0.0000327, 0.0017920],
          [0.0008414, -0.0004323, -0.0003952, -0.0005261, 0.0897077, 0.0004990],
          [-0.0013925, -0.0007727, -0.0015802, -0.0000313, -0.0003715,
           0.1960360]]  # INS-7382, matrice 7382M-1-1 pin 19 a 24

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

    return M1,M2,M3,M4,zeros1,zeros2,zeros3,zeros4

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

def Each_platform(ana,longueur,ana_vide,longueur_vide) :
    M1, M2, M3, M4, zeros1, zeros2, zeros3, zeros4 = matrices()
    ind_stop = Named_pins(c3d_experimental)
    longueur, ana, longueur_vide, ana_vide, longueur_statique, ana_statique = Affichage(c3d_experimental,c3d_vide,ind_stop)
    # vecteur des 6 forces/moments pour chaque plateforme

    platform1 = ana[0:6,:]
    platform2 = ana[6:12,:]
    platform3 = ana[12:18,:]
    platform4 = ana[18:24,:]

    platform1_vide, platform2_vide, platform3_vide, platform4_vide = plateforme_vide(c3d_vide)

    platform1_statique = ana_statique[0:6, :]
    platform2_statique = ana_statique[6:12, :]
    platform3_statique = ana_statique[12:18, :]
    platform4_statique = ana_statique[18:24, :]



    #donnees brutes statiques auxquelles on soustrait les valeurs a vide
    for j in range (6) :
        platform1_statique_brut[j,:]=platform1_statique[j,:] - np.mean(ana_vide1[j])
        platform2_statique_brut[j,:] = platform2_statique[j,:] - np.mean(ana_vide2[j])
        platform3_statique_brut[j,:] = platform3_statique[j,:] - np.mean(ana_vide3[j])
        platform4_statique_brut[j,:] = platform4_statique[j,:] - np.mean(ana_vide4[j])

    #donnees statiques auxquelles on a soustrait les valeurs a vide, puis calibre :
    platform1_statique_calib = np.matmul(M1, platform1_statique) * 100 * 10
    platform2_statique_calib = np.matmul(M2, platform2_statique) * 200 * 10
    platform3_statique_calib = np.matmul(M3, platform3_statique) * 100 * 10
    platform4_statique_calib = np.matmul(M4, platform4_statique) * 25 * 10




    return platform1_calib,platform2_calib,platform3_calib,platform4_calib




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
    plt.plot(x, (Poids_disque[numero_disque] / 4)*np.ones(longueur), '-g', label='Poids theorique divise par 4 pour ' + str(numero_disque) + ' disques')
    plt.legend()

if action == 'dynamique' :
    platform = dynamique(c3d_experimental)
    longueur = np.size(platform[0, 0])
    x = np.linspace(0, longueur, longueur)

    plt.subplot(311)
    plt.plot(x, platform[0, 0, :], '-k', label='Fx calibre PF1')
    plt.plot(x, platform[1, 0, :], '-b', label='Fx calibre PF2')
    plt.plot(x, platform[2, 0, :], '-r', label='Fx calibre PF3')
    plt.plot(x, platform[3, 0, :], '-m', label='Fx calibre PF4')
    plt.legend()
    plt.title('Essai dynamique ' + str(trial_name))

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
    plt.title('Fichier de calibration apres calibration : ' + str(calibrage_name_name))

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

plt.show()



