
import numpy as np
from ezc3d import c3d
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import seaborn as sns

# trial_name,date = 'statique_vide','dimanche'
# trial_name,date = 'statique_D8','dimanche'
# trial_name,date = 'statique_D11','dimanche'
# trial_name,date='SautHaut_03','dimanche'
trial_name,date='Test_Plateformes_4empilees05','samedi'
# trial_name,date = 'statiqueLeft_D8','dimanche'
# trial_name,date = 'statiqueLeft_D11','dimanche'
c3d_experimental = c3d(trial_name+'.c3d')

numero=2 #0:Fx, 1:Fy, 2:Fz, 3:Mx, 4:My, 5:Mz

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
    #on garde seulement les pins qui ont un nom
    force_labels=c3d_experimental['parameters']['ANALOG']['LABELS']['value']
    indices_supp=[]
    for i in range (len(force_labels)) :
        if 'Raw' in force_labels[i] :
            indices_supp=np.append(indices_supp,i)
    ind_stop=int(indices_supp[0])
    return ind_stop

def Affichage(c3d_experimental,ind_stop) :
    #affichage des donnees des plateformes
    if date=='samedi' :
        ana1=c3d_experimental['data']['analogs'][0,18:36,:] #F1x,F1y,F1z,M1x,M1y,M1z,.....M4z size 24
        ana2=c3d_experimental['data']['analogs'][0,45:51,:]
        ind=[i for i in range (18,36)]+[j for j in range (45,51)]
        ana=0
        longueur = np.size(ana1[0])
    else :
        ana=c3d_experimental['data']['analogs'][0,ind_stop:,:]
        longueur=np.size(ana[0])
        ana1,ana2=0,0
    # for i in range (4) :
    #     index=numero+6*i
    #     # longueur = np.size(ana[numero+6*i])
    #     print(str(c3d_experimental['parameters']['ANALOG']['LABELS']['value'][ind_stop:][index]) + ' = ' + str(
    #         ana[index, int(longueur / 2)]))
    return longueur,ana,ana1,ana2

def Each_platform(ana,ana1,ana2,longueur) :
    # vecteur des 6 forces/moments pour chaque plateforme
    zeros1 = np.array([1.0751899, 2.4828501, -0.1168980, 6.8177500, -3.0313399, -0.9456340])
    zeros2 = np.array([0., -2., -2., 0., 0., 0.])
    zeros3 = np.array([0.0307411, -5., -4., -0.0093422, -0.0079338, 0.0058189])
    zeros4 = np.array([-0.1032560, -3., -3., 0.2141770, 0.5169040, -0.3714130])


    platform1 = np.zeros((6, longueur))
    platform2 = np.zeros((6, longueur))
    platform3 = np.zeros((6, longueur))
    platform4 = np.zeros((6, longueur))

    if date=='samedi' :
        for i in range (longueur) :
            platform1[:,i] = ana1[0:6,i] + zeros1
            platform2[:,i] = ana1[6:12,i] + zeros2
            platform3[:,i] = ana1[12:18,i] + zeros3
            platform4[:,i] = ana2[0:6,i] + zeros4
    else :
        for i in range (longueur) :
            platform1[:,i] = ana[0:6,i] + zeros1
            platform2[:,i] = ana[6:12,i] + zeros2
            platform3[:,i] = ana[12:18,i] + zeros3
            platform4[:,i] = ana[18:24,i] + zeros4


    # M4 = [[-1.5038965, 0.0191295, -0.0143262, 0.0000000, -0.8255000, -0.0715982],
    #       [-0.0202159, -1.5082034, 0.0001771, 0.8255000, -0.0000000, 0.0582084],
    #       [-0.0132051, -0.0257038, -6.0646319, -0.0000000, -0.0000000, -0.0619103],
    #       [-0.0001964, 0.0002213, -0.0000134, 0.0100000, -0.0000000, -0.0014188],
    #       [-0.0000209, -0.0001944, 0.0001566, -0.0000000, 0.0100000, 0.0032266],
    #       [-0.0115609, -0.0099352, -0.0000278, 0.0000000, 0.0000000,-0.3888576]]  # Old_S2M, matrice Vieille_PF2_dimanche, pin 1 a 6

    # M4 = [[-0.4209, -0.1135,0.0044,-0.0013,0.017,-0.0127],
    #       [0.7716,0.4754,0.0062,-0.0008,-0.0095,-0.0026],
    #       [0.1891,0.2030,0.5321,0.0025,-0.0024,-0.0004],
    #       [0.0005,82.5496,0.,1.,-0.,-0.],
    #       [-82.5490,-0.0005,0.,-0.,1.,-0.],
    #       [-11.5812,3.5361,0.03,0.045,0.1378,0.3372]]  # Old_S2M, matrice Vieille_PF2_dimanche, pin 1 a 6

    # M4 = [[-20.7765, -0.3008,0.009,0.0037,0.2618,-0.0188],
    #       [11.0275,6.1685,0.0104,0.0647,-0.1339,-0.0024],
    #       [2.0835,1.5248,0.2088,0.0185,-0.0253,-0.001],
    #       [0.,82.5500,0.,1.,-0.,-0.],
    #       [-82.5489,-0.,0.,-0.,1.,-0.],
    #       [-43.2203,48.9740,0.0076,0.5886,0.5180,0.0614]]  # Old_S2M, matrice Vieille_PF2_dimanche, pin 1 a 6
    # M4 = [[5.7170, 0.1458, 0.0159, 0., -0.0461, 0.2619],
    #       [-2.2226, 2.1662, 0.0548, 0.0007, 0.0307, 0.0584],
    #       [0.0157, 0.0065, 0.1653, 0.0001, -0.0002, 0.0042],
    #       [-0., 0.3302, -0., 0.004, 0., 0.],
    #       [0.3302, 0., 0., 0., 0.004, -0.],
    #       [-4.8070, 1.3227, 0.0234, 0.0028, 0.0476, -0.0851]]

    # M4=[[7.5944, 0.1739,0.0115, 0., -0.0681, 0.2085],
    #       [-6.3583, 1.8617, 0.0466, 0., 0.0829, 0.2324],
    #       [0.0209, 0.0033, 0.1653, 0., -0.0002, 0.0049],
    #       [-0., 0.3302, -0., 0.004, 0., 0.],
    #       [-0.3302, -0., -0., -0., 0.004, 0.],
    #       [-0.4949, 0.0560, -0.0002, 0.0001, 0.0055, 0.0139]]\

    M4=[[4.1390,0.4018,-0.0456,-0.0280,-0.1606,-0.9432],
          [-0.7248,2.7855,0.1181,0.0009,0.1933,0.9311],
          [0.0111,-0.0063,0.8371,-0.004,0.0014,0.0120],
          [-4.4537,3.1811,-0.2881,0.9530,-1.2086,-4.5486],
          [-0.1021,0.0049,-0.0019,-0.0139,0.2170,-0.0016],
          [-0.1899,-0.0652,0.0083,0.0113,0.0341,0.1495]]

    # M4[0] = np.multiply(M4[0],-8)
    # M4[1]=np.multiply(M4[1],3.5)
    # M4[2]=np.multiply(M4[2],1.4)


    # M4=[[-157.4463,-1.3897,-0.6051,-0.0193,0.0146,4.1982],
    #     [6.4211,-155.9483,-2.1111,0.0041,-0.0275,0.4459],
    #     [5.2802,0.6117,-608.0808,-0.0031,-0.0106,-0.1563],
    #     [0.,82.5500,-0,1,-0,0],
    #     [-82.5500,-0,-0,-0,1,0],
    #     [-2.1179,6.8423,-5.6775,-0.1568,0.2920,-37,9273]]

    M1 = [[1.5289740, -0.0055372, 0.0000385, -0.0073467, 0.0008197, -0.7711510],
          [0.0145983, 1.5235620, 0.0000019, -0.0216603, -0.0003118, 0.7447910],
          [0.0041178, 0.0013416, 0.3924790, -0.0051392, -0.0000419, -0.0002830],
          [0.0077071, 0.0930658, 0.0254123, 1.4258270, 0.0021123, 0.8327660],
          [0.0059914, 0.0018655, -0.0010224, 0.0029237, 0.0917925, 0.0435810],
          [0.0017551, 0.0050230, -0.0020866, 0.0109248, -0.0002057, 0.1652510]]  # New_S2M, matrice 7273M pin 10 a 15


    M2 = [[1.5378351,    0.0014677,   -0.0003717,    0.0005122,   -0.0002392,   -0.0012910],
          [-0.0046930,    1.5345052,    0.0033831,    0.0002683,   -0.0001478,   -0.0006520],
          [-0.0016725,   -0.0021274,    0.3923144,    0.0001368,    0.0000125,   -0.0001880],
          [0.0019166,   -0.0004030,    0.0012135,    0.0898106,   -0.0000327,    0.0017920],
          [0.0008414,   -0.0004323,   -0.0003952,   -0.0005261,    0.0897077,    0.0004990],
          [-0.0013925,   -0.0007727,   -0.0015802,   -0.0000313,   -0.0003715,    0.1960360]]  # INS-7382, matrice 7382M-1-1 pin 19 a 24

    M3 = [[1.5438977, -0.0073794, 0.0022438, 0.0000509, -0.0002489, -0.0003690],
          [0.0104010, 1.5370655, 0.0013308, -0.0001873, 0.0000820, 0.0011140],
          [-0.0022109, -0.0003121, 0.3901421, 0.0001339, -0.0000556, -0.0002320],
          [0.0010562, -0.0010583, 0.0019478, 0.0900776, 0.0000017, 0.0014160],
          [0.0009278, -0.0021898, -0.0002942, 0.0002198, 0.0900668, 0.0003340],
          [-0.0003099, 0.0000449, 0.0035554, -0.0005566, -0.0002628, 0.1962430]]  # INS-7383, matrice 7383M-1-5 pin 28 a 33

    # M1,M2,M3,M4=np.transpose(M1),np.transpose(M2),np.transpose(M3),np.transpose(M4)

    platform1_calib = np.matmul(M1, platform1) * 100 * 10
    platform2_calib = np.matmul(M2, platform2) * 200 * 10
    platform3_calib = np.matmul(M3, platform3) * 100 * 10
    # platform4_calib = np.matmul(np.linalg.inv(M4), platform4) * 25 * 10
    platform4_calib = np.matmul(M4, platform4) * 25 * 10

    if date=='samedi' or 'Saut' in trial_name :
        zero_variable_1 = np.zeros(6)
        zero_variable_2 = np.zeros(6)
        zero_variable_3 = np.zeros(6)
        zero_variable_4 = np.zeros(6)

        for j in range (6) :
            zero_variable_1[j] = sum(platform1_calib[j, :100]) /100
            zero_variable_2[j] = sum(platform2_calib[j, :100]) / 100
            zero_variable_3[j] = sum(platform3_calib[j, :100]) / 100
            zero_variable_4[j] = sum(platform4_calib[j, :100]) / 100

        for i in range (longueur) :
            platform1_calib[:,i] = platform1_calib[:,i] - zero_variable_1
            platform2_calib[:,i] = platform2_calib[:,i] - zero_variable_2
            platform3_calib[:, i] = platform3_calib[:, i] - zero_variable_3
            platform4_calib[:, i] = platform4_calib[:, i] - zero_variable_4

    return platform1_calib,platform2_calib,platform3_calib,platform4_calib

def Courbes(ana,ana1,ana2,index,longueur) :
    platform1_calib, platform2_calib, platform3_calib, platform4_calib = Each_platform(ana,ana1,ana2,longueur)
    x = np.linspace(0, longueur, longueur)
    fig=plt.figure()
    fig.add_subplot(313)

    plt.subplot(311)
    plt.plot(x, platform1_calib[0, :], '-k', label='Fx calibre PF1')
    plt.plot(x, platform2_calib[0, :], '-b', label='Fx calibre PF2')
    plt.plot(x, platform3_calib[0, :], '-r', label='Fx calibre PF3')
    plt.plot(x, platform4_calib[0, :], '-m', label='Fx calibre PF4')
    plt.legend()
    plt.title('Essai ' + str(trial_name))

    plt.subplot(312)
    plt.plot(x, platform1_calib[1, :], '-k', label='Fy calibre  PF1')
    plt.plot(x, platform2_calib[1, :], '-b', label='Fy calibre PF2')
    plt.plot(x, platform3_calib[1, :], '-r', label='Fy calibre PF3')
    plt.plot(x, platform4_calib[1, :], '-m', label='Fy calibre PF4')
    plt.legend()

    plt.subplot(313)
    plt.plot(x, platform1_calib[2, :], '-k', label='Fz calibre PF1')
    plt.plot(x, platform2_calib[2, :], '-b', label='Fz calibre PF2')
    plt.plot(x, platform3_calib[2, :], '-r', label='Fz calibre PF3')
    plt.plot(x, platform4_calib[2, :], '-m', label='Fz calibre PF4')

    if 'D8' in trial_name:
        plt.plot(x,[1599.03/4 for i in range(longueur)],'-g',label='force theorique')
    if 'vide' in trial_name:
        plt.plot(x, [0 for i in range(longueur)], '-g', label='force theorique')
    if 'D11' in trial_name:
        plt.plot(x,[2242.566/4 for i in range(longueur)],'-g',label='force theorique')
    plt.legend()

    plt.show()

# plt.figure()
# plt.plot(c3d_experimental['data']['analogs'][0,-6:,:].T)
# plt.show()

ind_stop=Named_pins(c3d_experimental)
longueur,ana,ana1,ana2=Affichage(c3d_experimental,ind_stop)
Courbes(ana,ana1,ana2,numero,longueur)



