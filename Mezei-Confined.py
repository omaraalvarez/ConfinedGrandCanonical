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

import math as m
import numpy as np
import os

from random import random, randrange
from statistics import mean, pstdev
from timeit import default_timer as timer
           
def Mezei(ChemPot, L, h, T, R_Cut = 3.0):
    start = timer()
    """     AMOUNT OF STEPS     """
    MC_Relaxation_Steps = 2000
    MC_Equilibrium_Steps = 1000
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps
    MC_Measurement = 5

    x, y, z = [], [], []
    V = m.pow(L, 3.)
    z_Bins = 200
    Density_z = np.zeros(z_Bins)
    Beta = 1. / T
    Overlap = 0.75
    Pc, Pc_Sum, Pc_N = dict(), dict(), dict()
    Displacement, N_Displacement, N_Accepted_Displacement = L / 8, 0, 0
    N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
    N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
    N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0
    N_Measurements = 0
    Energy, Virial = 0., 0.
    Energy_Array, Pressure_Array, Density_Array = [], [], []

    Output_Route = "Output/ChPot_%.3f_h%.1f_T_%.2f" % (ChemPot, h, T)
    if not os.path.exists(Output_Route):
        os.makedirs(Output_Route)
    Average_Energy_File = open("%s/Average_Energy.dat" % Output_Route, "w+")
    Average_Energy_File.write("#\t< E / N >\n")
    Average_Pressure_File = open("%s/Average_Pressure.dat" % Output_Route, "w+")
    Average_Pressure_File.write("#\t< P >\n")
    Average_Density_File = open("%s/Average_Density.dat" % Output_Route, "w+")
    Average_Density_File.write("#\t< Density >\n")
    
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
            print("     Accepted: %d (%.2f)%%" % (N_Movement_Accepted, 100*N_Movement_Accepted/N_Movement))
            print("     Rejected: %d (%.2f)%%" % (N_Movement_Rejected, 100*N_Movement_Rejected/N_Movement))
            print("Insertions: %d" % N_Insertion)
            print("     Accepted: %d (%.2f)%%" % (N_Insertion_Accepted, 100*N_Insertion_Accepted/N_Insertion))
            print("     Rejected: %d (%.2f)%%" % (N_Insertion_Rejected, 100*N_Insertion_Rejected/N_Insertion))
            print("Removal: %d" % N_Removal)
            print("     Accepted: %d (%.2f)%%" % (N_Removal_Accepted, 100*N_Removal_Accepted/N_Removal))
            print("     Rejected: %d (%.2f)%%" % (N_Removal_Rejected, 100*N_Removal_Rejected/N_Removal))
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
            print("     Accepted: %d (%.2f)%%" % (N_Movement_Accepted, 100*N_Movement_Accepted/N_Movement))
            print("     Rejected: %d (%.2f)%%" % (N_Movement_Rejected, 100*N_Movement_Rejected/N_Movement))
            print("Insertions: %d" % N_Insertion)
            print("     Accepted: %d (%.2f)%%" % (N_Insertion_Accepted, 100*N_Insertion_Accepted/N_Insertion))
            print("     Rejected: %d (%.2f)%%" % (N_Insertion_Rejected, 100*N_Insertion_Rejected/N_Insertion))
            print("Removal: %d" % N_Removal)
            print("     Accepted: %d (%.2f)%%" % (N_Removal_Accepted, 100*N_Removal_Accepted/N_Removal))
            print("     Rejected: %d (%.2f)%%" % (N_Removal_Rejected, 100*N_Removal_Rejected/N_Removal))
            print("")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0

        if i == MC_Relaxation_Steps:
            print("~~~  STARTING MEASUREMENT STEPS  ~~~")

        if len(x) > 1:
            """            MOVEMENT           """
            N_Movement += 1
            N_Displacement += 1
            j = randrange(0, len(x), 1)
            Energy_Old, Virial_Old = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
            x_Old, y_Old, z_Old = x[j], y[j], z[j]

            """     DISPLACEMENT AND PERIODIC BOUNDARY CONDITIONS """
            x[j] += Displacement * (random() - 0.5)
            if x[j] <  -L / 2.:
                x[j] += L
            elif x[j] > L / 2.:
                x[j] -= L
                
            y[j] += Displacement * (random() - 0.5)
            if y[j] <  -L / 2.:
                y[j] += L
            elif y[j] > L / 2.:
                y[j] -= L

            z_Displacement = Displacement * (random() - 0.5)
            z[j] += z_Displacement
            if z[j] < - (h - 1) / 2.:
                z[j] -= z_Displacement
            elif z[j] > (h - 1) / 2.:
                z[j] -= z_Displacement
            
            Energy_New, Virial_New = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
            Delta_E = Energy_New - Energy_Old
            Delta_Virial = Virial_New - Virial_Old
            if random() < m.exp(-Beta * Delta_E):
                N_Movement_Accepted += 1
                N_Accepted_Displacement += 1
                Energy += Delta_E
                Virial += Delta_Virial
            else:
                N_Movement_Rejected += 1
                x[j] = x_Old
                y[j] = y_Old
                z[j] = z_Old

        RN = random()
        if RN < 0.5:
            """            INSERTION          """
            N_Insertion += 1
            grid = 8
            EVMPS = np.zeros((grid, grid, grid))
            Step = L / (grid + 1)
            Step_z = (h - 1) / (grid + 1)
            for i_x in range(0, grid, 1):
                for i_y in range(0, grid, 1):
                    for i_z in range(0, grid, 1):
                        x_grid = (i_x + 1) * Step - L / 2.
                        y_grid = (i_y + 1) * Step - L / 2.
                        z_grid = (i_z + 1) * Step_z - (h - 1) / 2.
                        for k in range(0, len(x), 1):
                            Delta_x = x_grid - x[k]
                            Delta_y = y_grid - y[k]
                            Delta_z = z_grid - z[k]
                            r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
                            if r2 < m.pow(Overlap, 2):
                                EVMPS[i_x, i_y, i_z] = 1
                                break
            if np.sum(EVMPS) < m.pow(grid, 3):
                aux1 = np.where(EVMPS == 0)
                j = randrange(0, len(aux1[0]), 1)
                i_x = aux1[0][j]
                i_y = aux1[1][j]
                i_z = aux1[2][j]
                x_Insertion = (i_x + 1) * Step - L / 2.
                y_Insertion = (i_y + 1) * Step - L / 2.
                z_Insertion = (i_z + 1) * Step_z - (h - 1) / 2.
                Energy_Insertion, Virial_Insertion = Energy_Virial(L, R_Cut, x_Insertion, y_Insertion, z_Insertion, x, y, z)
                """     PROBABILITY Pc      """
                #N_i = N_t
                #for j in range(N_t):
                #    rx = L*(random() - 0.5)
                #    ry = L*(random() - 0.5)
                #    rz = L*(random() - 0.5)
                #    for k in range(len(x)):
                #        Delta_x = rx - x[k]
                #        Delta_y = ry - y[k]
                #        Delta_z = rz - z[k]
                #        r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
                #        if r2 < m.pow(Overlap, 2):
                #            N_i -= 1
                            #if N_i < 1:
                            #    raise ValueError("Ni cant be negative")
                #            break

                if len(x) not in Pc_N:
                    Pc_N[len(x)] = 1
                    Pc_Sum[len(x)] = 1 - (np.sum(EVMPS) / m.pow(grid, 3)) #N_i/N_t
                    Pc[len(x)] = 1 - (np.sum(EVMPS) / m.pow(grid, 3)) #N_i/N_t
                else:
                    Pc_N[len(x)] += 1
                    Pc_Sum[len(x)] += 1 - (np.sum(EVMPS) / m.pow(grid, 3)) #N_i/N_t
                    Pc[len(x)] = Pc_Sum[len(x)] / Pc_N[len(x)]
                if random() < (V * Pc[len(x)] / (len(x) + 1)) * m.exp(Beta * (ChemPot - Energy_Insertion)):
                    N_Insertion_Accepted += 1
                    x.append(x_Insertion)
                    y.append(y_Insertion)
                    z.append(z_Insertion)
                    Energy += Energy_Insertion
                    Virial += Virial_Insertion
                else:
                    N_Insertion_Rejected += 1

                #"""     EXCLUDED VOLUME     """
                #V_Excluded = len(x) * (4. / 3.) * m.pi
                #""" CORRECTION DUE TO SPHERE OVERLAPNESS """
                #V_Excluded_Correction = 0
                #for j in range(0, len(x), 1):
                #    for k in range(0, len(x), 1):
                #        if j != k:
                #            Delta_x = x[j] - x[k]
                #            Delta_y = y[j] - y[k]
                #            Delta_z = z[j] - z[k]
                #            r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
                #            if r2 < 1.0: # r2 < 1 implies overlapness 
                #                d = m.sqrt(r2)
                #                V_Excluded_Correction += (m.pi / 12.) * (4 + d) * m.pow((2 - d), 2)
                #            else:
                #                if Delta_x > L - 1.: """  PERIODIC BOUNDARY CONDITION FOR THE OVERLAPNESS CONTRIBUTION  """
                #                    Delta_x -= L
                #                if Delta_y > L - 1.:
                #                    Delta_y -= L
                #                if Delta_z > L - 1.:
                #                    Delta_z -= L
                #                r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
                #                if r2 < 1.0: # r2 < 1 implies overlapness 
                #                    d = m.sqrt(r2)
                #                    V_Excluded_Correction += (m.pi / 12.) * (4 + d) * m.pow((2 - d), 2)
                #Volume_Ratio = (V_Excluded - V_Excluded_Correction) / V
                #if Volume_Ratio > 1 or Volume_Ratio < 0:
                #    raise ValueError("Volume Ratio (", Volume_Ratio, ") can't be negative nor greater than one.")
                #if random() < (V_Excluded - V_Excluded_Correction)/(len(x) + 1) * m.exp(Beta * (ChemPot - Energy_Insertion))
                #    N_Insertion_Accepted += 1
                #    x.append(x_Insertion)
                #    y.append(y_Insertion)
                #    z.append(z_Insertion)
                #    Energy += Energy_Insertion
                #else:
                #    N_Insertion_Rejected += 1

            else:        
                x_Insertion = L*(random() - 0.5)
                y_Insertion = L*(random() - 0.5)
                z_Insertion = (h - 1)*(random() - 0.5)
                Energy_Insertion, Virial_Insertion = Energy_Virial(L, R_Cut, x_Insertion, y_Insertion, z_Insertion, x, y, z)
                if random() < m.exp(Beta * (ChemPot - Energy_Insertion) + m.log(V / (len(x) + 1))):
                    N_Insertion_Accepted += 1
                    x.append(x_Insertion)
                    y.append(y_Insertion)
                    z.append(z_Insertion)
                    Energy += Energy_Insertion
                    Virial += Virial_Insertion
                else:
                    N_Insertion_Rejected += 1

        elif len(x) > 1:
            """            REMOVAL            """
            N_Removal += 1
            j = randrange(0, len(x), 1)

            """   IN CASE THE VALUE REQUIRED DOES'T EXIST
                    A LINEAR INTERPOLATION IS CARRIED ON.  """
            if (len(x) - 1) not in Pc:
                if len(Pc) == 1:
                    Pc_Interpolation =  Pc[ list( Pc.keys() )[0] ]
                elif len(Pc) > 1:

                    lim_inf = len(x) - 1
                    while True:
                        if lim_inf in Pc:
                            break
                        lim_inf -= 1

                    lim_sup = len(x) + 1
                    while True:
                        if lim_sup in Pc:
                            break
                        lim_sup += 1
                elif len(Pc) == 0:
                    raise ValueError("Pc has no values still.")

                Pc_Interpolation = Pc[lim_inf] * (1 - (len(x) - lim_inf) / (lim_sup - lim_inf)) + Pc[lim_sup] * ((len(x) - lim_inf) / (lim_sup - lim_inf))
            else:
                Pc_Interpolation = Pc[len(x) - 1]

            if random() > m.pow(1 - m.pow(Pc_Interpolation, len(x) - 1), m.pow(grid, 3)):
                Energy_Removal, Virial_Removal = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
                if random() < (len(x) / (V * Pc_Interpolation)) * m.exp(Beta * (Energy_Removal - ChemPot)):
                    N_Removal_Accepted += 1
                    x.pop(j)
                    y.pop(j)
                    z.pop(j)
                    Energy -= Energy_Removal
                    Virial -= Virial_Removal
                else:
                    N_Removal_Rejected += 1
            else: 
                Energy_Removal, Virial_Removal = Energy_Virial(L, R_Cut, x[j], y[j], z[j], x, y, z)
                if random() < m.exp(Beta * (Energy_Removal - ChemPot) + m.log(len(x) / V)):
                    N_Removal_Accepted += 1
                    x.pop(j)
                    y.pop(j)
                    z.pop(j)
                    Energy -= Energy_Removal
                    Virial -= Virial_Removal
                else:
                    N_Removal_Rejected += 1

        if i > 0:
            if m.fmod(i, MC_Measurement) == 0:
                if i >= MC_Relaxation_Steps:
                    """     COMPUTES AVERAGES    """
                    Energy_Array.append( Energy/len(x) )
                    Average_Energy_File.write("%d\t%f\n" % (len(Energy_Array), mean(Energy_Array)))
                    Pressure_Array.append( (len(x) * T - Virial / 3.0) / V )
                    Average_Pressure_File.write("%d\t%f\n" % (len(Pressure_Array), mean(Pressure_Array)))
                    Density_Array.append( len(x) / V ) 
                    Average_Density_File.write("%d\t%f\n" % (len(Density_Array), mean(Density_Array) ))
                    Density_z = z_Distribution(z, z_Bins, h, Density_z)
                    N_Measurements += 1
                """ PARTICLE DISPLACEMENT """
                if 1. * N_Accepted_Displacement / N_Displacement > 0.55:
                    Displacement *= 1.05
                else:
                    Displacement *= 0.95
                if Displacement < 0.05:
                    Displacement = 0.05
                if Displacement > L / 4.:
                    Displacement = L /4. 

                N_Accepted_Displacement, N_Displacement = 0, 0
    
    Average_Energy_File.close()
    Average_Pressure_File.close()
    Average_Density_File.close()

    Density_z = [(z_Bins * x) / (N_Measurements * V) for x in Density_z]
    Density_z_File = open("%s/Density_z.dat" % Output_Route, "w+")
    for i in range(len(Density_z)):
        Density_z_File.write("%.6f\t%.6f\n" % ((h - 1) * (i / z_Bins + (1 / z_Bins - 1) / 2), Density_z[i]))
    Density_z_File.close()

    Pc_File = open("%s/Pc.dat" % Output_Route, "w+")
    for i in Pc:
        Pc_File.write("%d\t%f\n" % (i, Pc[i]))
    Pc_File.close()

    print("< E / N > = %.6f + %.6f" % (mean(Energy_Array), pstdev(Energy_Array)) )
    print("< P > = %.6f + %.6f" % (mean(Pressure_Array), pstdev(Pressure_Array)) )
    print("< Density > = %.6f + %.6f " % (mean(Density_Array), pstdev(Density_Array)) )
    print("< N > = %.6f + %.6f " % (V*mean(Density_Array), V*pstdev(Density_Array)) )

    dt = timer() - start
    print("Elapsed time: %f s" % dt)

def u(r2):
    if r2 <= 1:
        return 100000
    if r2 > 1 and r2 <= 1.5:
        return -1
    if r2 > 1.5:
        return 0

#def u(r2):
#    return 4.0*(m.pow(1. / r2, 6.0) - m.pow(1. / r2, 3.))

def rdu(r2):
    return 4.0*(6*m.pow(1. / r2, 3.0) - 12.0 * m.pow(1.0 / r2, 6.0))

def Energy_Virial(L, R_Cut, rx, ry, rz, x, y, z):
    Energy, Virial = 0., 0.
    for i in range(len(x)):
        Delta_x = rx - x[i]
        Delta_y = ry - y[i]
        Delta_z = rz - z[i]        
        """    PERIODIC BOUNDARY CONDITIONS    """
        if Delta_x > L / 2.0:
            Delta_x -= L
        elif Delta_x < -L / 2.0:
            Delta_x += L

        if Delta_y > L / 2.0:
                Delta_y -= L
        elif Delta_y < -L / 2.0:
            Delta_y += L

        r2 = m.pow(Delta_x, 2) + m.pow(Delta_y, 2) + m.pow(Delta_z, 2)
        if r2 != 0.0:
            if r2 < m.pow(R_Cut, 2):
                Energy += u(r2)
                Virial += rdu(r2)
    return Energy, Virial

def z_Distribution(z, z_Bins, h, Density_z):
    z_sorted = z.copy()
    z_sorted.sort()
    Delta = (h - 1) / z_Bins
    for i in range(z_Bins):
        A = [x >= i * Delta - (h - 1) / 2 for x in z_sorted]
        B = [x <= (i + 1) * Delta - (h - 1) / 2 for x in z_sorted]
        Density_z[i] += np.sum(np.logical_and(A, B))
    return Density_z

Mezei(4.436, 10., 4., 2.0)
