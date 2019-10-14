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
    MC_Relaxation_Steps = 100000
    MC_Equilibrium_Steps = 5000000
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps
    MC_Measurement = 100
    """     VARIABLE    INITIALIZATION         """
    x, y, z = [], [], []
    V = h * pow(L, 2)
    Beta = 1. / T
    Pc_Random, Pc_Random_Sum, Pc_Random_N = dict(), dict(), dict()
    Displacement, N_Displacement, N_Displacement_Accepted = L / 8, 0, 0
    N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
    N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
    N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0
    Energy = 0.
    Energy_Array, Density_Array = [], []
    Bins = 100
    g_x = np.zeros(Bins, dtype = float)
    g_y = np.zeros(Bins, dtype = float)
    g_z = np.zeros(Bins, dtype = float)
    N_Measurements = 0

    """         PROGRAM     HEADER          """
    print("#" * 40 + "\n#" + " " * 38 + "#")
    print("#" + " " * 4 + "Monte Carlo + Excluded Volume" + " " * 5 + "#" + "\n#" + " " * 38 + "#")
    print("#" * 40)

    """       OUTPUT      ROUTES       """
    Output_Route = "Output/ChPot_%.3f_h_%.2f_L_%.2f_T_%.2f" % (ChemPot, h, L, T)
    if not os.path.exists(Output_Route):
        os.makedirs(Output_Route)
    if not os.path.exists("%s/Positions" % Output_Route):
        os.makedirs("%s/Positions" % Output_Route)
    Average_Energy_File = open("%s/Average_Energy.dat" % Output_Route, "w+")
    Average_Energy_File.write("#\t< E / N >\n")
    #Average_Pressure_File = open("%s/Average_Pressure.dat" % Output_Route, "w+")
    #Average_Pressure_File.write("#\t< P >\n")
    Average_Density_File = open("%s/Average_Density.dat" % Output_Route, "w+")
    Average_Density_File.write("#\t< Density >\n")

    """     SIMULATION'S    CYCLE       """
    for i in range(MC_Steps):    
        """ PRINTS TO SCREEN SUMMARY """
        if i > 0 and i < MC_Relaxation_Steps and i % int(.01 * MC_Relaxation_Steps) == 0:
            print("%d%% Relaxation" % (100*i / MC_Relaxation_Steps))
            print("U / N = %.6f" % (Energy/len(x)) )
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
            print("%d%% Equilibrium (%d Measurements)" % (100*(i - MC_Relaxation_Steps) / (MC_Equilibrium_Steps), N_Measurements))
            print("U / N = %.6f" % (Energy/len(x)) )
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
            Pc_File = open("%s/Pc_Relaxation.dat" % Output_Route, "w+")
            Pc_File.write("#N\t#Pc\n" )
            for i in Pc_Random:
                Pc_File.write("%d\t%.6f\n" % (i, Pc_Random[i]))
            Pc_File.close()

        RN = randrange(1, 4)
        """     EQUAL PROBABILITY FOR PARTICLE MOVEMENT, INSERTION OR DELETION      """
        if RN == 1 and len(x) > 1:
            """            MOVEMENT           """
            N_Movement += 1
            N_Displacement += 1
            Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted = Movement(h, L, Beta, Displacement, Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z)

        if RN == 2:
            """            INSERTION          """
            N_Insertion += 1
            Volume_Ratio, x_Insertion, y_Insertion, z_Insertion = Random_Excluded_Volume(h, L, x, y, z)
            if len(x) not in Pc_Random_N:
                Pc_Random_N[len(x)] = 1
                Pc_Random_Sum[len(x)] = Volume_Ratio
                Pc_Random[len(x)] = Volume_Ratio
            else:
                Pc_Random_N[len(x)] += 1
                Pc_Random_Sum[len(x)] += Volume_Ratio
                Pc_Random[len(x)] = Pc_Random_Sum[len(x)] / Pc_Random_N[len(x)]
            if len(x) > 0:
                """     INSERTION IN AVAILABLE SPACE    """
                Energy, N_Insertion_Accepted, N_Insertion_Rejected = Insertion_Mezei(h, L, V, Beta, ChemPot, R_Cut, Pc_Random, Pc_Random_Sum, Pc_Random_N, Energy, N_Insertion_Accepted, N_Insertion_Rejected, x, y, z, x_Insertion, y_Insertion, z_Insertion)                
            else:        
                """     RANDOM INSERTION        """
                Energy, N_Insertion_Accepted, N_Insertion_Rejected = Insertion(h, L, V, ChemPot, Beta, R_Cut, x, y, z, Energy, N_Insertion_Accepted, N_Insertion_Rejected)
        
        if RN == 3 and len(x) > 1:
            """            REMOVAL            """
            N_Removal += 1
            """   Pc[N - 1] IS REQUIRED  """
            if len(Pc_Random) == 1:
                Pc_Interpolation = list(Pc_Random.values())[0]
            else:
                if (len(x) - 1) not in Pc_Random:
                    """     IF Pc[N - 1] IS NOT IN THE ARRAY, AN EXTRAPOLATION OF ITS VALUE IS USED IN THE MEANTIME VIA THE CLOSEST HIGHER AND LOWER VALUE"""
                    Pc_Interpolation = Interpolation(Pc_Random, len(x))
                else:
                    Pc_Interpolation = Pc_Random[len(x) - 1]

            if random() > pow(1 - Pc_Interpolation, 1000):
                Energy, N_Removal_Accepted, N_Removal_Rejected = Removal_Mezei(Pc_Interpolation, h, L, V, Beta, ChemPot, R_Cut, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z)
            else: 
                Energy, N_Removal_Accepted, N_Removal_Rejected = Removal(h, L, V, Beta, ChemPot, R_Cut, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z)
        """      AVERAGES    """
        if i > 0:
            if fmod(i, MC_Measurement) == 0:
                if i >= MC_Relaxation_Steps:
                    Energy_Array.append( Energy/len(x) )
                    Average_Energy_File.write("%d\t%f\n" % (len(Energy_Array), mean(Energy_Array)))
                    Density_Array.append( len(x) / V ) 
                    Average_Density_File.write("%d\t%f\n" % (len(Density_Array), mean(Density_Array) ))
                    #g_x = Distribution(x, Bins, L, g_x)
                    #g_y = Distribution(y, Bins, L, g_y)
                    #g_z = Distribution(z, Bins, h, g_z)
                    g_x += Distribution(Bins, L, mean(Density_Array), x)
                    g_y += Distribution(Bins, L, mean(Density_Array), y)
                    g_z += Distribution(Bins, h, mean(Density_Array), z)
                    N_Measurements += 1
                    Positions_File = open("%s/Positions/Pos_%d.xyz" % (Output_Route, N_Measurements), "w+")
                    for i in range(len(x) - 1):
                        Positions_File.write("%.6f,\t%.6f,\t%.6f,\n" % (x[i], y[i], z[i]))
                    Positions_File.write("%.6f,\t%.6f,\t%.6f\n" % (x[i], y[i], z[i]))
                    Positions_File.close()
                """ PARTICLE DISPLACEMENT """
                if fmod(i, 100*MC_Measurement) == 0:
                    if 1. * N_Displacement_Accepted / N_Displacement > 0.55:
                        Displacement *= 1.05
                    else:
                        Displacement *= 0.95
                    if Displacement < 0.05:
                        Displacement = 0.05
                    if Displacement > L / 4.:
                        Displacement = L /4. 
                    N_Displacement, N_Displacement_Accepted = 0, 0
    
    Average_Energy_File.close()
    Average_Density_File.close()

    """     CAVITY PROBABILITIES OUTPUT     """
    Pc_File = open("%s/Pc.dat" % Output_Route, "w+")
    Pc_File.write("#N\t#Pc\n" )
    for i in Pc_Random:
        Pc_File.write("%d\t%.6f\n" % (i, Pc_Random[i]))
    Pc_File.close()

    """     Z DISTRIBUTION OUTPUT      """
    g_x = [Bins * x / (N_Measurements * V) for x in g_x]
    Delta = L / Bins
    g_x_File = open("%s/Density_x.dat" % Output_Route, "w+")
    for i in range(Bins):
        g_x_File.write( "%.6f\t%.6f\n" % (- L / 2 + (i + 0.5) * Delta, g_x[Bins - (i + 1)]) )
    g_x_File.close()

    g_y = [Bins * x / (N_Measurements * V) for x in g_y]
    Delta = L / Bins
    g_y_File = open("%s/Density_y.dat" % Output_Route, "w+")
    for i in range(Bins):
        g_y_File.write( "%.6f\t%.6f\n" % (- L / 2 + (i + 0.5) * Delta, g_y[Bins - (i + 1)]) )
    g_y_File.close()

    g_z = [Bins * x / (N_Measurements * V) for x in g_z]
    Delta = h / Bins
    g_z_File = open("%s/Density_z.dat" % Output_Route, "w+")
    for i in range(Bins):
        g_z_File.write( "%.6f\t%.6f\n" % (- h / 2 + (i + 0.5) * Delta, g_z[Bins - (i + 1)]) )
    g_z_File.close()

    print("< E / N > = %.6f + %.6f" % (mean(Energy_Array), pstdev(Energy_Array)) )
    print("< Density > = %.6f + %.6f " % (mean(Density_Array), pstdev(Density_Array)) )
    print("< N > = %.6f + %.6f " % (V*mean(Density_Array), V*pstdev(Density_Array)) )

    Summary_File = open("%s/Summary.txt" % Output_Route, "w+")
    Summary_File.write("Mezei algorithm for the Grand Canonical Monte Carlo Simulation\n\n")
    Summary_File.write("~~~~~   INPUT   ~~~~~\n\n")
    Summary_File.write("Chemical Potential = %.6f\n" % ChemPot)
    Summary_File.write("Volume = %.3f (L = %.6f)\n" % (V, L))
    Summary_File.write("Temperature = %.3f\n\n" % T)
    Summary_File.write("Relaxation Steps: %d.   Equilibrium Steps: %d\n\n" % (MC_Relaxation_Steps, MC_Equilibrium_Steps))
    Summary_File.write("~~~~~   OUTPUT  ~~~~~\n\n")
    Summary_File.write("< E / N > = %.6f + %.6f\n" % (mean(Energy_Array), pstdev(Energy_Array)) )
    Summary_File.write("< Density > = %.6f + %.6f\n" % (mean(Density_Array), pstdev(Density_Array)) )
    Summary_File.write("< N > = %.6f + %.6f\n" % (V*mean(Density_Array), V*pstdev(Density_Array)) )

    dt = timer() - start
    print("Elapsed time: %f s" % dt)

def Random_Excluded_Volume(h, L, x, y, z):
    N_Random = 200
    N_in = 0
    x_Insertion, y_Insertion, z_Insertion = [], [], []
    if not len(x) > 1:
        for _ in range(N_Random):
            z_V = h*(random() - 0.5)
            if h / 2 - abs(z_V) > 0.5:
                N_in += 1
                x_V = L*(random() - 0.5)
                y_V = L*(random() - 0.5)
                x_Insertion.append(x_V)
                y_Insertion.append(y_V)
                z_Insertion.append(z_V)
    else:
        for _ in range(N_Random):
            z_V = h*(random() - 0.5)
            if h / 2 - abs(z_V) > 0.5:
                for j in range(len(x)):
                    x_V = L*(random() - 0.5)
                    y_V = L*(random() - 0.5)
                    Delta_x = x_V - x[j]
                    Delta_x = PeriodicBoundaryConditions(L, Delta_x)
                    Delta_y = y_V - y[j]
                    Delta_y = PeriodicBoundaryConditions(L, Delta_y)
                    Delta_z = z_V - z[j]
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
    
def Movement(h, L, Beta, Displacement, Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z):
    j = randrange(0, len(x), 1)
    Energy_Old = Energy_Calculation(h, L, R_Cut, x[j], y[j], z[j], x, y, z)
    x_Old, y_Old, z_Old = x[j], y[j], z[j]
    x[j] += Displacement * (random() - 0.5)                
    x[j] = PeriodicBoundaryConditions(L, x[j])
    y[j] += Displacement * (random() - 0.5)
    y[j] = PeriodicBoundaryConditions(L, y[j])
    z[j] += Displacement * (random() - 0.5)
    Energy_New = Energy_Calculation(h, L, R_Cut, x[j], y[j], z[j], x, y, z)
    Delta_E = Energy_New - Energy_Old
    if random() < exp(-Beta * Delta_E):
        N_Movement_Accepted += 1
        N_Displacement_Accepted += 1
        Energy += Delta_E
    else:
        N_Movement_Rejected += 1
        x[j] = x_Old
        y[j] = y_Old
        z[j] = z_Old
    return Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted

def Insertion_Mezei(h, L, V, Beta, ChemPot, R_Cut, Pc, Pc_Sum, Pc_N, Energy, N_Insertion_Accepted, N_Insertion_Rejected, x, y, z, x_Insertion, y_Insertion, z_Insertion):
    j = randrange(0, len(x_Insertion), 1)
    Energy_Insertion = Energy_Calculation(h, L, R_Cut, x_Insertion[j], y_Insertion[j], z_Insertion[j], x, y, z)
    if random() < (V * Pc[len(x)] / (len(x) + 1)) * exp(Beta * (ChemPot - Energy_Insertion)):
        N_Insertion_Accepted += 1
        x.append(x_Insertion[j])
        y.append(y_Insertion[j])
        z.append(z_Insertion[j])
        Energy += Energy_Insertion
    else:
        N_Insertion_Rejected += 1
    return Energy, N_Insertion_Accepted, N_Insertion_Rejected

def Insertion(h, L, V, ChemPot, Beta, R_Cut, x, y, z, Energy, N_Insertion_Accepted, N_Insertion_Rejected):
    x_Insertion = L*(random() - 0.5)
    y_Insertion = L*(random() - 0.5)
    z_Insertion = h*(random() - 0.5)
    Energy_Insertion = Energy_Calculation(h, L, R_Cut, x_Insertion, y_Insertion, z_Insertion, x, y, z)
    if random() < exp(Beta * (ChemPot - Energy_Insertion) + log(V / (len(x) + 1))):
        N_Insertion_Accepted += 1
        x.append(x_Insertion)
        y.append(y_Insertion)
        z.append(z_Insertion)
        Energy += Energy_Insertion
    else:
        N_Insertion_Rejected += 1
    return Energy, N_Insertion_Accepted, N_Insertion_Rejected

def Removal_Mezei(Pc_Interpolation, h, L, V, Beta, ChemPot, R_Cut, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z):
    j = randrange(0, len(x), 1)
    Energy_Removal = Energy_Calculation(h, L, R_Cut, x[j], y[j], z[j], x, y, z)
    if random() < (len(x) / (V * Pc_Interpolation)) * exp(Beta * (Energy_Removal - ChemPot)):
        N_Removal_Accepted += 1
        x.pop(j)
        y.pop(j)
        z.pop(j)
        Energy -= Energy_Removal
    else:
        N_Removal_Rejected += 1
    return Energy, N_Removal_Accepted, N_Removal_Rejected

def Removal(h, L, V, Beta, ChemPot, R_Cut, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z):
    j = randrange(0, len(x), 1)
    Energy_Removal = Energy_Calculation(h, L, R_Cut, x[j], y[j], z[j], x, y, z)
    if random() < exp(Beta * (Energy_Removal - ChemPot) + log(len(x) / V)):
        N_Removal_Accepted += 1
        x.pop(j)
        y.pop(j)
        z.pop(j)
        Energy -= Energy_Removal
    else:
        N_Removal_Rejected += 1
    return Energy, N_Removal_Accepted, N_Removal_Rejected

def u_ParticleParticle(r2):
    """     SQUARE WELL POTENTIAL       """
    if r2 <= 1:
        return inf
    elif r2 <= pow(1.5, 2):
        return -1
    else:
        return 0
    """     LENNARD JONES POTENTIAL     """
    #return 4.0*(pow(1. / r2, 6.0) - pow(1. / r2, 3.))

def u_ParticleWall(h, z):
    Delta_z = h / 2 - abs(z)
    """     SQUARE WELL POTENTIAL       """
    if Delta_z <= 0.5:
        return inf
    elif Delta_z <= 1.5:#h / 2:
        #if z > 0:
        return -1
        #else:
        #    return 1
    else:
        return 0 

def PeriodicBoundaryConditions(L, x):
    if x <  -L / 2.:
        x += L
    elif x > L / 2.:
        x -= L
    
    return x

def Energy_Calculation(h, L, R_Cut, rx, ry, rz, x, y, z):
    Energy = u_ParticleWall(h, rz)
    if Energy == inf:
        return Energy 
    for i in range(len(x)):
        Delta_x = rx - x[i]
        Delta_y = ry - y[i]
        Delta_z = rz - z[i]        
        Delta_x = PeriodicBoundaryConditions(L, Delta_x)
        Delta_y = PeriodicBoundaryConditions(L, Delta_y)
        r2 = pow(Delta_x, 2) + pow(Delta_y, 2) + pow(Delta_z, 2)
        if r2 != 0.0:
            if r2 < pow(R_Cut, 2):
                Energy += u_ParticleParticle(r2)
    return Energy

def Interpolation(Pc, l):
    """     LINEAR INTERPOLATION BY TAKING THE PREVIOUS AND FOLLOWING VALUE """
    if len(Pc) == 1:
        Pc_Interpolation =  Pc[ list( Pc.keys() )[0] ]
    elif len(Pc) > 1:
        lim_inf = l - 1
        while True:
            if lim_inf in Pc:
                break
            lim_inf -= 1

        lim_sup = l + 1
        while True:
            if lim_sup in Pc:
                break
            lim_sup += 1
    elif len(Pc) == 0:
        raise ValueError("Pc has no values still.")
    Pc_Interpolation = Pc[lim_inf] * ((1 - l - lim_inf) / (lim_sup - lim_inf)) + Pc[lim_sup] * ((l - lim_inf) / (lim_sup - lim_inf))
    return Pc_Interpolation

def Distribution(Bins, h, Density, z):
    Delta = h / Bins
    g_z = np.zeros(Bins, dtype = float)
    for i in range( len(z) ): 
        l = int( (h / 2 - z[i]) / Delta) - 1
        g_z[l] += 1
    #for l in range(Bins):
    #    g_z[l] = g_z[l] / ( len(z) )
    return g_z

Mezei(-3.0, 10.0, 5., 2.0)