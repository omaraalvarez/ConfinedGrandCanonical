##############################################################
#                                                            #
#       Mezei's algorithm for dense fluids from              #
#       Molecular Physics: An International Journal          #
#       at the Interface Between Chemistry and               #
#       Physics: A cavity-biased (T, V, μ) Monte Carlo       #
#       method for the computer simulation of                #
#       fluids. (1980)                                       #
#                                                            #
##############################################################

import numpy as np
import os

from math import exp, log, pow, fmod, pi, sqrt, ceil, inf
from random import random, randrange
from statistics import mean, pstdev
from timeit import default_timer as timer

def Mezei(ChemPot, h, L, T, R_Cut = 3.0):
    start = timer()
    """     CONFIGURATIONAL     STEPS           """
    MC_Relaxation_Steps = 15000
    MC_Equilibrium_Steps = 25000
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps
    MC_Measurement = 10
    """     VARIABLE    INITIALIZATION         """
    x, y, z = [], [], []
    V = h * pow(L, 2)
    Beta = 1. / T
    Pc_Random, Pc_Random_Sum, Pc_Random_N = dict(), dict(), dict()
    Displacement, N_Displacement, N_Displacement_Accepted = L / 8, 0, 0
    N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
    N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
    N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0
    Energy, Virial = 0., 0.
    Energy_Array, Pressure_Array, Density_Array = [], [], []
    z_Bins = 200
    N_Measurements = 0

    """         PROGRAM     HEADER          """
    print("#" * 40 + "\n#" + " " * 38 + "#")
    print("#" + " " * 4 + "Monte Carlo + Excluded Volume" + " " * 5 + "#" + "\n#" + " " * 38 + "#")
    print("#" * 40)

        """       OUTPUT      ROUTES       """
    Output_Route = "Output/ChPot_%.3f_h_%.2f_T_%.2f" % (ChemPot, h, T)
    if not os.path.exists(Output_Route):
        os.makedirs(Output_Route)
    Average_Energy_File = open("%s/Average_Energy.dat" % Output_Route, "w+")
    Average_Energy_File.write("#\t< E / N >\n")
    Average_Pressure_File = open("%s/Average_Pressure.dat" % Output_Route, "w+")
    Average_Pressure_File.write("#\t< P >\n")
    Average_Density_File = open("%s/Average_Density.dat" % Output_Route, "w+")
    Average_Density_File.write("#\t< Density >\n")

    """     SIMULATION'S    CYCLE       """
    for i in range(MC_Steps):    
        """ PRINTS TO SCREEN SUMMARY """
        if i > 0 and i < MC_Relaxation_Steps and i % int(.01 * MC_Relaxation_Steps) == 0:
            print("%d%% Relaxation" % (100*i / MC_Relaxation_Steps))
            print("U / N = %.6f" % (Energy/len(x)) )
            print("P = %.6f" % ((len(x) * T  - Virial / 3.0) / V) )
            print("N = %d" % len(x))
            print("Density = %.6f" % (len(x)/V) )
            print("Max Displacement = %.6f" % Displacement)
            print("Movements: %d" % N_Movement)
            print("     Accepted: %d" % N_Movement_Accepted)
            print("     Rejected: %d" % N_Movement_Rejected)
            print("Insertions: %d" % N_Insertion)
            print("     Accepted: %d" % N_Insertion_Accepted)
            print("     Rejected: %d" % N_Insertion_Rejected)
            print("Removal: %d" % N_Removal)
            print("     Accepted: %d" % N_Removal_Accepted)
            print("     Rejected: %d" % N_Removal_Rejected)
            print("")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0

        if i > 0 and i > MC_Relaxation_Steps and i % int(0.01 * (MC_Steps - MC_Relaxation_Steps)) == 0:
            print("%d%% Equilibrium" % (100*(i - MC_Relaxation_Steps) / (MC_Equilibrium_Steps)))
            print("U / N = %.6f" % (Energy/len(x)) )
            print("P = %.6f" % ((len(x) * T  - Virial / 3.0) / V) )
            print("N = %d" % len(x))
            print("Density = %.6f" % (len(x)/V) )
            print("Max Displacement = %.6f" % Displacement)
            print("Movements: %d" % N_Movement)
            print("     Accepted: %d" % N_Movement_Accepted)
            print("     Rejected: %d" % N_Movement_Rejected)
            print("Insertions: %d" % N_Insertion)
            print("     Accepted: %d" % N_Insertion_Accepted)
            print("     Rejected: %d" % N_Insertion_Rejected)
            print("Removal: %d" % N_Removal)
            print("     Accepted: %d" % N_Removal_Accepted)
            print("     Rejected: %d" % N_Removal_Rejected)
            print("")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0

        if i == MC_Relaxation_Steps:
            print("~~~  STARTING MEASUREMENT STEPS  ~~~")

        RN = randrange(1, 4)
        """     EQUAL PROBABILITY FOR PARTICLE MOVEMENT, INSERTION OR DELETION      """

        if RN == 2:
            """            INSERTION          """
            N_Insertion += 1
            Volume_Ratio, x_Insertion, y_Insertion, z_Insertion = Random_Excluded_Volume(h, L, x, y, z)



def Random_Excluded_Volume(h, L, x, y, z):
    N_Random = 1000
    N_in = 0
    x_Insertion, y_Insertion, z_Insertion = [], [], []
    if not len(x) > 1:
        for _ in range(N_Random):
            x_V = L*(random() - 0.5)
            y_V = L*(random() - 0.5)
            z_V = h*(random() - 0.5)
            N_in += 1
            x_Insertion.append(x_V)
            y_Insertion.append(y_V)
            z_Insertion.append(z_V)
    else:
        for _ in range(N_Random):
            x_V = L*(random() - 0.5)
            y_V = L*(random() - 0.5)
            z_V = h*(random() - 0.5)
            for j in range(len(x)):
                Delta_x = x_V - x[j]
                Delta_x = PeriodicBoundaryConditions(L, Delta_x)
                Delta_y = y_V - y[j]
                Delta_y = PeriodicBoundaryConditions(L, Delta_y)
                Delta_z = z_V - z[j]
                Delta_z = PeriodicBoundaryConditions(L, Delta_z)
                r2 = pow(Delta_x, 2) + pow(Delta_y, 2) + pow(Delta_z, 2)
                if r2 < pow(0.7, 2):
                    break
                if j == len(x) - 1:
                    N_in += 1
                    x_Insertion.append(x_V)
                    y_Insertion.append(y_V)
                    z_Insertion.append(z_V)
    Volume_Ratio = N_in / N_Random
    return Volume_Ratio, x_Insertion, y_Insertion, z_Insertion