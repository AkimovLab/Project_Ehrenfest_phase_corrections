import matplotlib
matplotlib.use('Agg')

import sys
import matplotlib.pyplot as plt
#from matplotlib.pyplot import figure
import numpy as np

"""
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
"""

#plt.rc('axes', titlesize=18)      # fontsize of the axes title
#plt.rc('axes', labelsize=14)      # fontsize of the x and y labels
#plt.rc('legend', fontsize=18)    # legend fontsize

plt.rc('axes', titlesize=24)      # fontsize of the axes title
plt.rc('axes', labelsize=20)      # fontsize of the x and y labels
plt.rc('legend', fontsize=20)     # legend fontsize
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels

plt.rc('figure.subplot', left=0.2)
plt.rc('figure.subplot', right=0.95)
plt.rc('figure.subplot', bottom=0.13)
plt.rc('figure.subplot', top=0.88)




#plt.rc('font', size=24)          # controls default text sizes
#plt.rc('figure', titlesize=24)    # fontsize of the figure title


print(plt.style.available)
#plt.style.use('ggplot')
#plt.style.use('fivethirtyeight')
#plt.style.use('dark_background')

#sys.exit(0)

colors = {}

colors.update({"11": "#8b1a0e"})  # red       
colors.update({"12": "#FF4500"})  # orangered 
colors.update({"13": "#B22222"})  # firebrick 
colors.update({"14": "#DC143C"})  # crimson   

colors.update({"21": "#5e9c36"})  # green
colors.update({"22": "#006400"})  # darkgreen  
colors.update({"23": "#228B22"})  # forestgreen
colors.update({"24": "#808000"})  # olive      

colors.update({"31": "#8A2BE2"})  # blueviolet
colors.update({"32": "#00008B"})  # darkblue  

colors.update({"41": "#2F4F4F"})  # darkslategray




def plot_2D_pes(model):


    prefix = "2D/model%i/" % (model)

    q, h0, h1, e0, e1 = np.loadtxt(prefix+"_pes.txt", usecols=(0, 1, 2, 3, 4), unpack=True)

    plt.title('Diabatic PES')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Energy, a.u.')
    plt.plot(q, h0, label='$H_{00}$', linewidth=4, color = colors["11"]) 
    plt.plot(q, h1, label='$H_{11}$', linewidth=4, color = colors["21"]) 
    plt.legend()
    plt.show()
    figname = prefix+"_pes-dia.png"
    plt.savefig(figname, dpi=300)
    plt.close()



    plt.title('Adiabatic PES')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Energy, a.u.')
    plt.plot(q, e0, label='$E_{0}$', linewidth=4, color = colors["11"]) 
    plt.plot(q, e1, label='$E_{1}$', linewidth=4, color = colors["21"]) 
    plt.legend()
    plt.show()
    figname = prefix+"_pes-adi.png"
    plt.savefig(figname, dpi=300)
    plt.close()



def plot_3D_pes(model):


    prefix = "3D/model%i/" % (model)

    q, h0, h1, h2, e0, e1, e2 = np.loadtxt(prefix+"_pes.txt", usecols=(0, 1, 2, 3, 4, 5, 6), unpack=True)

    plt.title('Diabatic PES')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Energy, a.u.')
    plt.plot(q, h0, label='$H_{00}$', linewidth=4, color = colors["11"]) 
    plt.plot(q, h1, label='$H_{11}$', linewidth=4, color = colors["21"]) 
    plt.plot(q, h2, label='$H_{22}$', linewidth=4, color = colors["32"]) 
    plt.legend()
    plt.show()
    figname = prefix+"_pes-dia.png"
    plt.savefig(figname, dpi=300)
    plt.close()



    plt.title('Adiabatic PES')
    plt.xlabel('Coordinate, a.u.')
    plt.ylabel('Energy, a.u.')
    plt.plot(q, e0, label='$E_{0}$', linewidth=4, color = colors["11"]) 
    plt.plot(q, e1, label='$E_{1}$', linewidth=4, color = colors["21"]) 
    plt.plot(q, e2, label='$E_{2}$', linewidth=4, color = colors["32"]) 
    plt.legend()
    plt.show()
    figname = prefix+"_pes-adi.png"
    plt.savefig(figname, dpi=300)
    plt.close()



def plot_2D(model, case, reord, phase):

    # X - P - Y
    # X - is the representation in which the dynamics is run
    # Y - is the population of a given state in that representation
    prefix = "2D/model%i/case%i_do_reordering%i_do_phase_correction%i_rep_" % (model, case, reord, phase)

    t0, dia_P_adi, dia_P_dia = np.loadtxt(prefix+"0_new.txt", usecols=(0, 9, 11), unpack=True)
    t1, adi_P_adi, adi_P_dia = np.loadtxt(prefix+"1_new.txt", usecols=(0, 9, 11), unpack=True)

    dia_q, dia_p, dia_Ekin, dia_Epot, dia_Etot = np.loadtxt(prefix+"0_new.txt", usecols=(1, 2,  3, 4, 5), unpack=True)
    adi_q, adi_p, adi_Ekin, adi_Epot, adi_Etot = np.loadtxt(prefix+"1_new.txt", usecols=(1, 2,  3, 4, 5), unpack=True)


    #================== Diabatic populaiton ============================
    plt.title('Diabatic Population of P(0)')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population')
    plt.plot(t0, dia_P_dia, label='P(0), Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_P_dia, label='P(0), Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "2D/model%i/_fig-%i-%i-%i-pop_dia.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Adiabatic populaiton ============================
    plt.title('Adiabatic Population of P(0)')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population')
    plt.plot(t0, dia_P_adi, label='P(0), Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_P_adi, label='P(0), Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "2D/model%i/_fig-%i-%i-%i-pop_adi.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Coordinate ============================
    plt.title('Coordinate')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Position, a.u.')
    plt.plot(t0, dia_q, label='Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_q, label='Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "2D/model%i/_fig-%i-%i-%i-t-q.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Momentum ============================
    plt.title('Momentum')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Momentum, a.u.')
    plt.plot(t0, dia_p, label='Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_p, label='Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "2D/model%i/_fig-%i-%i-%i-t-p.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Phase space diagram ============================
    plt.title('Phase space portrait')
    plt.xlabel('Position, a.u.')
    plt.ylabel('Momentum, a.u.')
    plt.plot(dia_q, dia_p, label='Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(adi_q, adi_p, label='Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "2D/model%i/_fig-%i-%i-%i-q-p.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Total energy ============================
    plt.title('Total energy')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Energy, a.u.')
    plt.plot(t0, dia_Etot, label='Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_Etot, label='Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "2D/model%i/_fig-%i-%i-%i-t-Etot.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()




def plot_3D(model, case, reord, phase):

    # X - P - Y
    # X - is the representation in which the dynamics is run
    # Y - is the population of a given state in that representation
    prefix = "3D/model%i/case%i_do_reordering%i_do_phase_correction%i_rep_" % (model, case, reord, phase)

    t0, dia_P_adi, dia_P_dia = np.loadtxt(prefix+"0_new.txt", usecols=(0, 9, 11), unpack=True)
    t1, adi_P_adi, adi_P_dia = np.loadtxt(prefix+"1_new.txt", usecols=(0, 9, 11), unpack=True)

    dia_q, dia_p, dia_Ekin, dia_Epot, dia_Etot = np.loadtxt(prefix+"0_new.txt", usecols=(1, 2,  3, 4, 5), unpack=True)
    adi_q, adi_p, adi_Ekin, adi_Epot, adi_Etot = np.loadtxt(prefix+"1_new.txt", usecols=(1, 2,  3, 4, 5), unpack=True)

    dia_cum_phase = np.loadtxt(prefix+"0_new.txt", usecols=(13), unpack=True)
    adi_cum_phase = np.loadtxt(prefix+"1_new.txt", usecols=(13), unpack=True)

    dia_phase = np.loadtxt(prefix+"0_new.txt", usecols=(17), unpack=True)
    adi_phase = np.loadtxt(prefix+"1_new.txt", usecols=(17), unpack=True)



    #================== Diabatic populaiton ============================
    plt.title('Diabatic Population of P(0)')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population')
    plt.plot(t0, dia_P_dia, label='P(0), Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_P_dia, label='P(0), Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "3D/model%i/_fig-%i-%i-%i-pop_dia.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Adiabatic populaiton ============================
    plt.title('Adiabatic Population of P(0)')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Population')
    plt.plot(t0, dia_P_adi, label='P(0), Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_P_adi, label='P(0), Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "3D/model%i/_fig-%i-%i-%i-pop_adi.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Coordinate ============================
    plt.title('Coordinate')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Position, a.u.')
    plt.plot(t0, dia_q, label='Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_q, label='Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "3D/model%i/_fig-%i-%i-%i-t-q.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Momentum ============================
    plt.title('Momentum')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Momentum, a.u.')
    plt.plot(t0, dia_p, label='Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_p, label='Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "3D/model%i/_fig-%i-%i-%i-t-p.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Phase space diagram ============================
    plt.title('Phase space portrait')
    plt.xlabel('Position, a.u.')
    plt.ylabel('Momentum, a.u.')
    plt.plot(dia_q, dia_p, label='Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(adi_q, adi_p, label='Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "3D/model%i/_fig-%i-%i-%i-q-p.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Total energy ============================
    plt.title('Total energy')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Energy, a.u.')
    plt.plot(t0, dia_Etot, label='Dyn = dia', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_Etot, label='Dyn = adi', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "3D/model%i/_fig-%i-%i-%i-t-Etot.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()

    #================== Cumulative phase ============================
    plt.title('Phases')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Phase')
    plt.plot(t0, adi_cum_phase, label='Cum. phase', linewidth=4, color = colors["21"]) 
    plt.plot(t1, adi_phase, label='Phase corr.', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "3D/model%i/_fig-%i-%i-%i-t-phase.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()


    plt.title('Cumulative phase')
    plt.xlabel('Time, a.u.')
    plt.ylabel('Phase')
    plt.plot(t0, adi_cum_phase, label='With correction', linewidth=4, color = colors["21"]) 
    plt.plot(t1, dia_cum_phase, label='No correction', linewidth=1, color = colors["11"]) 
    plt.legend()
    plt.show()
    figname = "3D/model%i/_fig-%i-%i-%i-t-phase-compare.png" % (model, case, reord, phase)
    plt.savefig(figname, dpi=300)
    plt.close()






def plot_toc():

    # X - P - Y
    # X - is the representation in which the dynamics is run
    # Y - is the population of a given state in that representation
    model = 1
    case = 0
    reord = 0
    phase = 0

    prefix = "2D/model%i/case%i_do_reordering%i_do_phase_correction%i_rep_" % (model, case, reord, phase)
    t0, dia_P_adi, dia_P_dia = np.loadtxt(prefix+"0_new.txt", usecols=(0, 9, 11), unpack=True)
    t1, adi_P_adi, adi_P_dia = np.loadtxt(prefix+"1_new.txt", usecols=(0, 9, 11), unpack=True)

    dia_q, dia_p, dia_Ekin, dia_Epot, dia_Etot = np.loadtxt(prefix+"0_new.txt", usecols=(1, 2,  3, 4, 5), unpack=True)
    adi_q, adi_p, adi_Ekin, adi_Epot, adi_Etot = np.loadtxt(prefix+"1_new.txt", usecols=(1, 2,  3, 4, 5), unpack=True)


    prefix = "2D/model%i/" % (model)
    q, h0, h1, e0, e1 = np.loadtxt(prefix+"_pes.txt", usecols=(0, 1, 2, 3, 4), unpack=True)




    #================== Adiabatic populaiton ============================
 

    plt.rc('axes', titlesize=10)     # fontsize of the axes title
    plt.rc('axes', labelsize=4)      # fontsize of the x and y labels
    plt.rc('legend', fontsize=9)    # legend fontsize
    plt.rc('font', size=0)          # controls default text sizes
    plt.rc('xtick', labelsize=4)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=4)    # fontsize of the tick labels
    plt.rc('figure', titlesize=0)    # fontsize of the figure title



    fig = plt.figure(num=None, figsize=(2, 2), dpi=300, frameon=False)

    ax1 = fig.add_subplot(111)
    ax1.plot(t0, dia_P_adi, label='Consistent Phases', linewidth=2, color = colors["21"]) 
    ax1.plot(t1, adi_P_adi, label='Inconsistent Phases', linewidth=2, color = colors["11"]) 
    ax1.tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False, labelbottom=False) 


    plt.title('Population Dynamics')
#    plt.xlabel('Time')
#    plt.ylabel('Population')
    ax1.legend(loc=8, bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True) # facecolor="#9d9e8e", edgecolor="b")
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.85, wspace=0.2, hspace=0.0)


    plt.show()
    plt.savefig("TOC.png", dpi=300)
    plt.close()
#    plt.clf()







for model in [1, 2, 3, 4, 5]:

    plot_2D_pes(model)

    for case in [0, 1, 2]:
        for do_reordering in [0, 1]:
            for do_phase_correction in [0, 1]:    
#                pass
                plot_2D(model, case, do_reordering, do_phase_correction)



for model in [1, 2, 3, 4, 5]:

    plot_3D_pes(model)

    for case in [0, 1, 2, 3]:
        for do_reordering in [0, 1]:
            for do_phase_correction in [0, 1]:    
#                pass
                plot_3D(model, case, do_reordering, do_phase_correction)

plot_toc()

