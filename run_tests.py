#*********************************************************************************
#* Copyright (C) 2017-2018 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

"""
  This file demonstrates how to run an ensemble of Ehrenfest trajectories
  using the hierarchy of Hamiltonians approach and the built in Ehrenfest1
  function
 
"""


import cmath
import math
import os
import sys


if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *






def compute_model(q, params, full_id):

    model = params["model"]
    res = None

    Id = Cpp2Py(full_id)
    indx = Id[-1]

    if model==1:
        res = models_Libra.model1(q.col(indx), params)
    elif model==2:
        res = models_Libra.model2(q.col(indx), params)
    elif model==3:
        res = models_Libra.model3(q.col(indx), params)
    elif model==4:
        res = models_Libra.model4(q.col(indx), params)

    elif model==6:
        res = models_Libra.model6(q.col(indx), params)
    elif model==7:
        res = models_Libra.model7(q.col(indx), params)


#    res.rep = params["rep"]    

    return res
    




def run_Ehrenfest(ndia, nadi, nnucl, ntraj, _q, _p, _Cdia, _Cadi, _iM, rep, outname, params):
    """
    model - setup the Hamiltonian
    rep - 0 - diabatic, 1 - adiabatic
    outname - the name of the output file
    """

    # Create copies of the input dynamical variables, so we could run several run_test 
    # functions with the same input variables without worries that they will be altered
    # inside of run_test

    q = MATRIX(_q)
    p = MATRIX(_p)
    iM = MATRIX(_iM)
    Cdia = CMATRIX(_Cdia)
    Cadi = CMATRIX(_Cadi)

    params.update( {"rep": rep } )

    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.init_all(2)
    print "id=", ham.id, " level=", ham.level


    ham1 = [] 
    for tr in xrange(ntraj):
        ham1.append( nHamiltonian(ndia, nadi, nnucl) )        
        print ham1[tr].id, ham1[tr].level
        ham1[tr].init_all(2)
        ham.add_child(ham1[tr])
        print Cpp2Py(ham1[tr].get_full_id())


    #  Set up the models and compute internal variables

    # Initial calculations
    ham.compute_diabatic(compute_model, q, params, 1)
    ham.compute_adiabatic(1, 1); 
    ham.ampl_adi2dia(Cdia, Cadi, 0, 1)

    Cdia.show_matrix()
    Cadi.show_matrix()

    if rep==0:
        ham.compute_nac_dia(p, iM, 0, 1);  
        ham.compute_hvib_dia(1); 
    elif rep==1:
        ham.compute_nac_adi(p, iM, 0, 1);
        ham.compute_hvib_adi(1); 


    ham1[0].get_nac_adi().show_matrix()
    ham1[0].get_ham_adi().show_matrix()
    ham1[0].get_ham_dia().show_matrix()

    Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh.compute_etot(ham, p, Cdia, Cadi, iM, rep)


#    cum_phase0 = ham1[0].get_cum_phase_corr().get(0)
#    cum_phase1 = ham1[0].get_cum_phase_corr().get(1)


    cum_phase0 = ham.get_cum_phase_corr(Py2Cpp_int([0,0])).get(0)
    cum_phase1 = ham.get_cum_phase_corr(Py2Cpp_int([0,0])).get(1)


    out = open(outname, "w")
    out.close()

    # Do the propagation
    dt = params["dt"]
    for i in xrange(params["nsteps"]):

        if rep==0:
            Ehrenfest1(dt, q, p, iM, Cdia, ham, compute_model, params, rep)

        elif rep==1:
#            Ehrenfest1(dt, q, p, iM, Cadi, ham, compute_model, params, rep)
            Ehrenfest2(dt, q, p, iM, Cadi, ham, compute_model, params, rep, params["do_reordering"], params["do_phase_correction"])



        #=========== Properties ==========
        dm_dia, dm_adi = tsh.compute_dm(ham, Cdia, Cadi, rep, 1)

        Ekin, Epot, Etot, dEkin, dEpot, dEtot = tsh.compute_etot(ham, p, Cdia, Cadi, iM, rep)


#        phase0 = ham1[0].get_cum_phase_corr().get(0) / cum_phase0
#        phase1 = ham1[0].get_cum_phase_corr().get(1) / cum_phase1

#        cum_phase0 = ham1[0].get_cum_phase_corr().get(0)
#        cum_phase1 = ham1[0].get_cum_phase_corr().get(1)

        phase0 = 1.0+0.0j
        if abs(cum_phase0) > 0.0:
            phase0 = ham.get_cum_phase_corr(Py2Cpp_int([0,0])).get(0) / cum_phase0

        phase1 = 1.0+0.0j
        if abs(cum_phase1) > 0.0:
            phase1 = ham.get_cum_phase_corr(Py2Cpp_int([0,0])).get(1) / cum_phase1

        cum_phase0 = ham.get_cum_phase_corr(Py2Cpp_int([0,0])).get(0)
        cum_phase1 = ham.get_cum_phase_corr(Py2Cpp_int([0,0])).get(1)


        out = open(outname, "a")
        ret = (i*dt, q.get(0), p.get(0), 
               Ekin, Epot, Etot, dEkin, dEpot, dEtot,
               dm_adi.get(0,0).real, dm_adi.get(1,1).real, dm_dia.get(0,0).real, dm_dia.get(1,1).real,
               cum_phase0.real, cum_phase0.imag,  cum_phase1.real, cum_phase1.imag,
               phase0.real, phase0.imag,   phase1.real, phase1.imag
              )
        out.write("%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f \
                   %8.5f %8.5f  %8.5f %8.5f    %8.5f %8.5f  %8.5f %8.5f\n" %  ret )
        out.close()




def run_test2(params, case, do_reordering, do_phase_correction): 

    ndia, nadi, nnucl, ntraj = 2, 2, 1, 1

    params.update( { "do_reordering": do_reordering, "do_phase_correction":do_phase_correction})

    rnd = Random()

    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   mean_q.set(0,0, params["q0"])
    sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, params["sq"])
    mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, params["p0"])
    sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, params["sp"])

    q = MATRIX(nnucl,ntraj);  tsh.sample(q, mean_q, sigma_q, rnd)
    p = MATRIX(nnucl,ntraj);  tsh.sample(p, mean_p, sigma_p, rnd)
    iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/params["mass"])

    Cdia, Cadi = CMATRIX(ndia, ntraj), CMATRIX(nadi, ntraj)
    for traj in xrange(ntraj):
        if case == 0:
            Cadi.set(0, traj, 1.0+0.0j);
        elif case == 1:
            Cadi.set(1, traj, 1.0+0.0j);
        elif case == 2:
            Cadi.set(0, traj, 1.0+0.0j);
            Cadi.set(1, traj, 1.0+0.0j);
            Cadi *= (1.0/math.sqrt(2.0))  


    prefix = params["prefix"] + "case"+str(case)+"_do_reordering"+str(do_reordering)+"_do_phase_correction"+str(do_phase_correction)+"_rep"

    run_Ehrenfest(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, 0, prefix+"_0_new.txt", params)
    run_Ehrenfest(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, 1, prefix+"_1_new.txt", params)



def run_test3(params, case, do_reordering, do_phase_correction): 

    ndia, nadi, nnucl, ntraj = 3, 3, 1, 1

    params.update( { "do_reordering": do_reordering, "do_phase_correction":do_phase_correction})

    rnd = Random()

    # Dynamical variables and system-specific properties
    mean_q = MATRIX(nnucl,1);   mean_q.set(0,0, params["q0"])
    sigma_q = MATRIX(nnucl,1);  sigma_q.set(0,0, params["sq"])
    mean_p = MATRIX(nnucl,1);   mean_p.set(0,0, params["p0"])
    sigma_p = MATRIX(nnucl,1);  sigma_p.set(0,0, params["sp"])

    q = MATRIX(nnucl,ntraj);  tsh.sample(q, mean_q, sigma_q, rnd)
    p = MATRIX(nnucl,ntraj);  tsh.sample(p, mean_p, sigma_p, rnd)
    iM = MATRIX(nnucl,1);     iM.set(0,0, 1.0/params["mass"])

    Cdia, Cadi = CMATRIX(ndia, ntraj), CMATRIX(nadi, ntraj)
    for traj in xrange(ntraj):
        if case == 0:
            Cadi.set(0, traj, 1.0+0.0j);
        elif case == 1:
            Cadi.set(1, traj, 1.0+0.0j);
        elif case == 2:
            Cadi.set(2, traj, 1.0+0.0j);
        elif case == 3:
            Cadi.set(0, traj, 1.0+0.0j);
            Cadi.set(1, traj, 1.0+0.0j);
            Cadi.set(2, traj, 1.0+0.0j);
            Cadi *= (1.0/math.sqrt(3.0))  


    prefix = params["prefix"] + "case"+str(case)+"_do_reordering"+str(do_reordering)+"_do_phase_correction"+str(do_phase_correction)+"_rep"

    run_Ehrenfest(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, 0, prefix+"_0_new.txt", params)
    run_Ehrenfest(ndia, nadi, nnucl, ntraj, q, p, Cdia, Cadi, iM, 1, prefix+"_1_new.txt", params)



def run2D(params):
    for do_reordering in [0, 1]:
        for do_phase_correction in [0, 1]:
            for case in [0, 1, 2]:
                print "=======", case, do_reordering, do_phase_correction, "=========="
                run_test2(params, case, do_reordering, do_phase_correction)


def run3D(params):
    for do_reordering in [0, 1]:
        for do_phase_correction in [0, 1]:
            for case in [0, 1, 2, 3]:
                print "=======", case, do_reordering, do_phase_correction, "=========="
                run_test3(params, case, do_reordering, do_phase_correction)



def model_profile(params, x0, dx, npts):

    nadi, ndia, nnucl, ntraj = 1, 1, 1, 1

    if params["model"]==6:  # 2D
        nadi, ndia = 2, 2 
    elif params["model"]==7:  # 3D
        nadi, ndia = 3, 3

    H, E = CMATRIX(ndia, ndia), CMATRIX(nadi, nadi)
   
    # ======= Hierarchy of Hamiltonians =======
    ham = nHamiltonian(ndia, nadi, nnucl)
    ham.init_all(2)

    ham1 = [] 
    for tr in xrange(ntraj):
        ham1.append( nHamiltonian(ndia, nadi, nnucl) )        
        print ham1[tr].id, ham1[tr].level
        ham1[tr].init_all(2)
        ham.add_child(ham1[tr])
        print Cpp2Py(ham1[tr].get_full_id())

    # Bind the Hamiltonian with the matrices:
    ham1[0].set_ham_dia_by_ref(H)
    ham1[0].set_ham_adi_by_ref(E)



    f = open(params["prefix"]+"_pes.txt", "w"); f.close()
    q = MATRIX(1,1);  

    for i in xrange(npts):

        q.set(0,0, x0 + dx*i)

        ham.compute_diabatic(compute_model, q, params, 1)
        ham.compute_adiabatic(1, 1); 

        f = open(params["prefix"]+"_pes.txt", "a")
        if params["model"]==6:  # 2D
            f.write( "%8.5f   %8.5f   %8.5f   %8.5f   %8.5f \n" % (q.get(0), H.get(0,0).real, H.get(1,1).real, E.get(0,0).real, E.get(1,1).real ))

        elif params["model"]==7:  # 3D
            f.write( "%8.5f   %8.5f   %8.5f   %8.5f  %8.5f   %8.5f   %8.5f\n" % 
            (q.get(0), H.get(0,0).real, H.get(1,1).real, H.get(2,2).real, E.get(0,0).real, E.get(1,1).real, E.get(2,2).real ))
        f.close()




x0, dx, npts = -5.0, 0.1, 1000


params = { "model":6 }
params.update( { "A0":0.1, "A1": 0.1 } )
params.update( { "B0":0.0, "B1": -0.1 } )
params.update( { "w0":0.25, "w1":0.25 } )
params.update( { "delta0":0.0, "delta1": 0.5*math.pi } )
params.update( { "V01": 0.01 } )
params.update( { "nsteps":10000, "dt":1.0 })
params.update( { "q0":0.1, "sq":0.00, "p0":0.0, "sp":0.00, "mass":2000.0 })


params.update( { "prefix": "2D/model1/"})
model_profile(params, x0, dx, npts)
run2D(params)

params.update( { "q0":0.1, "sq":0.00, "p0":5.0, "sp":0.00, "mass":2000.0 })
params.update( { "prefix": "2D/model2/"})
model_profile(params, x0, dx, npts)
run2D(params)

params.update( { "q0":0.1, "sq":0.00, "p0":25.0, "sp":0.00, "mass":2000.0 })
params.update( { "prefix": "2D/model3/"})
model_profile(params, x0, dx, npts)
run2D(params)

params.update( { "B1": 0.2 } )
params.update( { "q0":0.1, "sq":0.00, "p0":0.0, "sp":0.00, "mass":2000.0 })
params.update( { "prefix": "2D/model4/"})
model_profile(params, x0, dx, npts)
run2D(params)


params.update( { "q0":0.1, "sq":0.00, "p0":25.0, "sp":0.00, "mass":2000.0 })
params.update( { "prefix": "2D/model5/"})
model_profile(params, x0, dx, npts)
run2D(params)





params = { "model":7 }
params.update( { "A0":0.1, "A1": -0.1, "A2":-0.2} )
params.update( { "B0":0.0, "B1": 0.1, "B2":0.15 } )
params.update( { "w0":0.25, "w1":0.25, "w2":0.20 } )
params.update( { "delta0":0.0, "delta1": 0.0, "delta2":0.0 } )
params.update( { "V01": 0.05, "V02":0.05, "V12":0.03 } )
params.update( { "nsteps":10000, "dt":1.0})


params.update( { "q0":0.1, "sq":0.00, "p0":0.0, "sp":0.00, "mass":2000.0 })
params.update( { "prefix": "3D/model1/"})
model_profile(params, x0, dx, npts)
run3D(params)

params.update( { "q0":0.1, "sq":0.00, "p0":5.0, "sp":0.00, "mass":2000.0 })
params.update( { "prefix": "3D/model2/"})
model_profile(params, x0, dx, npts)
run3D(params)

params.update( { "q0":0.1, "sq":0.00, "p0":25.0, "sp":0.00, "mass":2000.0 })
params.update( { "prefix": "3D/model3/"})
model_profile(params, x0, dx, npts)
run3D(params)


params.update( { "q0":20.0, "sq":0.00, "p0":0.0, "sp":0.00, "mass":2000.0 })
params.update( { "prefix": "3D/model5/"})
model_profile(params, x0, dx, npts)
run3D(params)



params.update( { "A0":0.1, "A1": -0.1, "A2":-0.2} )
params.update( { "B0":0.0, "B1": 0.5, "B2":1.5 } )
params.update( { "w0":0.25, "w1":0.25, "w2":0.20 } )
params.update( { "delta0":0.0, "delta1": 0.0, "delta2":0.0 } )
params.update( { "V01": 0.05, "V02":0.05, "V12":0.03 } )
params.update( { "nsteps":10000, "dt":1.0})

params.update( { "q0":0.1, "sq":0.00, "p0":0.0, "sp":0.00, "mass":2000.0 })
params.update( { "prefix": "3D/model4/"})
model_profile(params, x0, dx, npts)
run3D(params)



