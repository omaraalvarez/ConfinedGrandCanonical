using Plots; gr();
using Plots.PlotMeasures;
using StatsPlots;
using Statistics;
using Distributions;
using Test;
using CSV;

function Pore_Separation()
    ChemPot = -4.150918;
    L = 20.;
    T = 2.;
    Bulk_Density = 0.1;
    H = range(2., 10., step = 1.);
    σ_p, λ_p = 0.5, 1.5;
    σ_w, λ_w = 0.5, 1.5;

    Patch_Percentage = 0;
    Patch_Radius = round(sqrt(L^2 * Patch_Percentage / 314) , digits = 6)

    Mean_Energy = zeros(Float64, length(H));
    Std_Energy = zeros(Float64, length(H));
    Mean_Density = zeros(Float64, length(H));
    Std_Density = zeros(Float64, length(H));
    j = 1;

    Output_Route = pwd() * "/Output_UPDATING/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/"
    mkpath("$Output_Route")
    p_x, p_y, p_z = Array{Any, 1}(undef, length(H)), Array{Any, 1}(undef, length(H)), Array{Any, 1}(undef, length(H));
    p_Normalized_x, p_Normalized_y, p_Normalized_z = Array{Any, 1}(undef, length(H)), Array{Any, 1}(undef, length(H)), Array{Any, 1}(undef, length(H));
    p_Energy, p_Histogram_Energy, p_Average_Energy = Array{Any, 1}(undef, length(H)), Array{Any, 1}(undef, length(H)), Array{Any, 1}(undef, length(H));
    p_Density, p_Histogram_Density, p_Average_Density = Array{Any, 1}(undef, length(H)), Array{Any, 1}(undef, length(H)), Array{Any, 1}(undef, length(H));

    TotalEnergy_File = open("$Output_Route/Energy_T_$(T)_ChemPot_$(round(ChemPot, digits = 2))_Patch_$Patch_Percentage%.dat", "w");
    TotalDensity_File = open("$Output_Route/Density_T_$(T)_ChemPot_$(round(ChemPot, digits = 2))_Patch_$Patch_Percentage%.dat", "w");
    println(TotalEnergy_File, "h\tEnergy\tStandardDeviation")
    println(TotalDensity_File, "h\tDensity\tStandardDeviation") 
    for h in H
        Mean_Energy[j], Std_Energy[j], Mean_Density[j], Std_Density[j], p_Energy[j], p_Histogram_Energy[j], p_Average_Energy[j], p_Density[j], p_Histogram_Density[j], p_Average_Density[j], p_x[j], p_y[j], p_z[j], p_Normalized_x[j], p_Normalized_y[j], p_Normalized_z[j] = MarkovChainMonteCarlo_Mezei(ChemPot, h, L, T, σ_p, λ_p, σ_w, λ_w, j, length(H), Patch_Percentage, Patch_Radius, Bulk_Density);
        println(TotalDensity_File, "$h\t$(Mean_Density[j])\t$(Std_Density[j])")
        println(TotalEnergy_File, "$h\t$(Mean_Energy[j])\t$(Std_Energy[j])")
        j += 1;
    end
    close(TotalEnergy_File)
    close(TotalDensity_File)

    Pore_Plots("X_Distribution", "X Axis", "Density", p_x[1], p_x[2], p_x[3], p_x[4], p_x[5], p_x[6], p_x[7], p_x[8], p_x[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Y_Distribution", "Y Axis", "Density", p_y[1], p_y[2], p_y[3], p_y[4], p_y[5], p_y[6], p_y[7], p_y[8], p_y[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Z_Distribution", "Z Axis", "Density", p_z[1], p_z[2], p_z[3], p_z[4], p_z[5], p_z[6], p_z[7], p_z[8], p_z[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Normalized_X_Distribution", "X Axis", "Normalized Density", p_Normalized_x[1], p_Normalized_x[2], p_Normalized_x[3], p_Normalized_x[4], p_Normalized_x[5], p_Normalized_x[6], p_Normalized_x[7], p_Normalized_x[8], p_Normalized_x[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Normalized_Y_Distribution", "Y Axis", "Normalized Density", p_Normalized_y[1], p_Normalized_y[2], p_Normalized_y[3], p_Normalized_y[4], p_Normalized_y[5], p_Normalized_y[6], p_Normalized_y[7], p_Normalized_y[8], p_Normalized_y[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Normalized_Z_Distribution", "Z Axis", "Normalized Density", p_Normalized_z[1], p_Normalized_z[2], p_Normalized_z[3], p_Normalized_z[4], p_Normalized_z[5], p_Normalized_z[6], p_Normalized_z[7], p_Normalized_z[8], p_Normalized_z[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Energy", "Measurements", "Energy [Unitless]", p_Energy[1], p_Energy[2], p_Energy[3], p_Energy[4], p_Energy[5], p_Energy[6], p_Energy[7], p_Energy[8], p_Energy[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Energy_Histogram", "Energy [Unitless]", "Normalized Frequency", p_Histogram_Energy[1], p_Histogram_Energy[2], p_Histogram_Energy[3], p_Histogram_Energy[4], p_Histogram_Energy[5], p_Histogram_Energy[6], p_Histogram_Energy[7], p_Histogram_Energy[8], p_Histogram_Energy[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Energy_Average", "Measurements", "Average Energy [Unitless]", p_Average_Energy[1], p_Average_Energy[2], p_Average_Energy[3], p_Average_Energy[4], p_Average_Energy[5], p_Average_Energy[6], p_Average_Energy[7], p_Average_Energy[8], p_Average_Energy[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Density", "Measurements", "Density", p_Density[1], p_Density[2], p_Density[3], p_Density[4], p_Density[5], p_Density[6], p_Density[7], p_Density[8], p_Density[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Density_Histogram", "Density", "Normalized Frequency", p_Histogram_Density[1], p_Histogram_Density[2], p_Histogram_Density[3], p_Histogram_Density[4], p_Histogram_Density[5], p_Histogram_Density[6], p_Histogram_Density[7], p_Histogram_Density[8], p_Histogram_Density[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    Pore_Plots("Density_Average", "Measurements", "Average Density", p_Average_Density[1], p_Average_Density[2], p_Average_Density[3], p_Average_Density[4], p_Average_Density[5], p_Average_Density[6], p_Average_Density[7], p_Average_Density[8], p_Average_Density[9], Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)


    N_Slab = 2 * convert(Int64, ceil(L / (2σ_w) + 3 + 1) * ceil(L / (√3 * 2σ_w) + √3))
    File_Slabs_Povray = open("$Output_Route/Slab.xyz", "w+")
    for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
        i_x == ceil(L / (2σ_w) + 3 + 1) && i_y == ceil(L / (√3 * 2σ_w) + √3) ? println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w),\t0.,\n$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w),\t0.") : println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w),\t0.,\n$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w),\t0.,")
    end
    close(File_Slabs_Povray)

    for h in H
        run(`povray $Output_Route/h_$h/Positions/Pore_X_Axis.ini`)
        run(`povray $Output_Route/h_$h/Positions/Pore_Z_Axis.ini`)
    end
    run(`povray $Output_Route/Pore_Slab.pov +W800 +H800`)

end

function MarkovChainMonteCarlo_Mezei(ChemPot::Float64, h::Float64, L::Float64, T::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, Number_Run::Int64, Total_Run::Int64, Patch_Percentage::Int64, Patch_Radius::Float64, Bulk_Density::Float64, R_Cut::Float64 = 3.)
    ####################################################     CONFIGURATIONAL STEPS   ################################################################
    V = (h + 2σ_w) * L^2;
    MC_Measurement = convert(Int64, ceil( 8V ));
    MC_Relaxation_Measurement = 1_000;
    MC_Equilibrium_Measurement = 10_000;
    MC_Relaxation_Steps = MC_Measurement * MC_Relaxation_Measurement;
    MC_Equilibrium_Steps = MC_Measurement * MC_Equilibrium_Measurement;
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps;

    ####################################################     ARRAY  INITIALIZATION   ################################################################
    N_Bins, Beta = 200, 1. / T;
    x, y, z = Float64[], Float64[], Float64[];
    Pc_Random, Pc_Random_Sum, Pc_Random_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    #Pc_Grid, Pc_Grid_Sum, Pc_Grid_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    #Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    Max_Displacement_Z, Max_Displacement, N_Displacement, N_Displacement_Accepted, N_Image = (h + 2σ_w) / 8., L / 8., 0, 0, 1;
    N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
    N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
    N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
    Energy, N_Measurements = 0., 0;
    Energy_Array, Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    Mean_Energy_Array, Mean_Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    Std_Energy, Std_Density = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ); 
    g_x, g_y, g_z = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ), N_Bins), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ), N_Bins), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ), N_Bins);
    PotentialFunction = zeros(Float64, N_Bins);

    ####################################################     OUTPUT ROUTE AND FILES   ################################################################
    Output_Route = pwd() * "/Output_UPDATING/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/h_$(h)"
    mkpath("$Output_Route/Positions")
    Acceptance_File = open("$Output_Route/Acceptance.dat", "w");
    println(Acceptance_File, "N\tMovement\tInsertion\tRemoval")
    ####################################################     SLAB CREATION   ################################################################
    File_Slabs_Povray = open("$Output_Route/Positions/Slab.xyz", "w+")
    for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
        if i_x == ceil(L / (2σ_w) + 3 + 1) && i_y == ceil(L / (√3 * 2σ_w) + √3)
            println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w),\t$(h / 2 + σ_w),")
            println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w),\t$(h / 2 + σ_w)")
        else
            println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w),\t$(h / 2 + σ_w),")
            println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w),\t$(h / 2 + σ_w),")
        end
    end
    close(File_Slabs_Povray)

    ####################################################    SIMULATION LOOP   ################################################################
    @inbounds for i = 1:MC_Steps
        ####################################################    RELAXATION STEPS SUMMARY   ################################################################
        if i < MC_Relaxation_Steps && i % floor(0.01MC_Relaxation_Steps) == 0
            @test all(Array(x) .<= L / 2.) && all(Array(x) .>= -L / 2.)
            @test all(Array(y) .<= L / 2.) && all(Array(y) .>= -L / 2.)
            @test all(Array(z) .<= h / 2.) && all(Array(z) .>= -h / 2.)
            println("$(convert(Int64, floor(100i / MC_Relaxation_Steps)))% Relaxation. [$Number_Run / $Total_Run]")
            println("μ = $ChemPot\tT = $T\tρ_Bulk = $Bulk_Density")
            println("Steps: $i / $MC_Relaxation_Steps")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("N = $(length(x))")
            println("Density = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Max_Displacement, digits = 6))")
            println("Max Displacement Z = $(round(Max_Displacement_Z, digits = 6))")
            println("Patch Radius: $Patch_Radius\tPatch Percetage: $Patch_Percentage%")
            println("Movements: $N_Movement")
            println("   Accepted: $N_Movement_Accepted ($(round(100N_Movement_Accepted / N_Movement, digits = 2))%)")
            println("   Rejected: $N_Movement_Rejected ($(round(100N_Movement_Rejected / N_Movement, digits = 2))%)")
            println("Insertions: $N_Insertion")
            println("   Accepted: $N_Insertion_Accepted ($(round(100N_Insertion_Accepted / N_Insertion, digits = 2))%)")
            println("   Rejected: $N_Insertion_Rejected ($(round(100N_Insertion_Rejected / N_Insertion, digits = 2))%)")
            println("Removal: $N_Removal")
            println("   Accepted: $N_Removal_Accepted ($(round(100N_Removal_Accepted / N_Removal, digits = 2))%)")
            println("   Rejected: $N_Removal_Rejected ($(round(100N_Removal_Rejected / N_Removal, digits = 2))%)")
            println("")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
        end

        ####################################################    EQUILIBRIUM STEPS SUMMARY   ################################################################
        if i > MC_Relaxation_Steps && i % floor(0.01MC_Equilibrium_Steps) == 0
            @test all(Array(x) .<= L / 2.) && all(Array(x) .>= -L / 2.)
            @test all(Array(y) .<= L / 2.) && all(Array(y) .>= -L / 2.)
            @test all(Array(z) .<= h / 2.) && all(Array(z) .>= -h / 2.)
            println("$(convert(Int64, floor(100(i - MC_Relaxation_Steps) / MC_Equilibrium_Steps)))% Equilibrium. [$Number_Run / $Total_Run]")
            println("μ = $ChemPot\tT = $T\tρ_Bulk = $Bulk_Density")
            println("Steps: $(i - MC_Relaxation_Steps) / $MC_Equilibrium_Steps.\t[$N_Measurements Measurements]")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("N = $(length(x))")
            println("Density = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Max_Displacement, digits = 6))")
            println("Max Displacement Z = $(round(Max_Displacement_Z, digits = 6))")
            println("Patch Radius: $Patch_Radius\tPatch Percetage: $Patch_Percentage%")
            println("Movements: $N_Movement")
            println("   Accepted: $N_Movement_Accepted ($(round(100N_Movement_Accepted / N_Movement, digits = 2))%)")
            println("   Rejected: $N_Movement_Rejected ($(round(100N_Movement_Rejected / N_Movement, digits = 2))%)")
            println("Insertions: $N_Insertion")
            println("   Accepted: $N_Insertion_Accepted ($(round(100N_Insertion_Accepted / N_Insertion, digits = 2))%)")
            println("   Rejected: $N_Insertion_Rejected ($(round(100N_Insertion_Rejected / N_Insertion, digits = 2))%)")
            println("Removal: $N_Removal")
            println("   Accepted: $N_Removal_Accepted ($(round(100N_Removal_Accepted / N_Removal, digits = 2))%)")
            println("   Rejected: $N_Removal_Rejected ($(round(100N_Removal_Rejected / N_Removal, digits = 2))%)")
            println("")
            println(Acceptance_File, "$(convert(Int64, floor(100(i - MC_Relaxation_Steps) / MC_Equilibrium_Steps)))\t$(round(100N_Movement_Accepted / N_Movement, digits = 2))\t$(round(100N_Insertion_Accepted / N_Insertion, digits = 2))\t$(round(100N_Removal_Accepted / N_Removal, digits = 2))")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
        end

        ####################################################    SAVES POSITIONS OF THE MOLECULES   ################################################################
        if i > MC_Relaxation_Steps && i % .1MC_Equilibrium_Steps == 0
            Positions_File = open("$Output_Route/Positions/Pos_$N_Image.xyz", "w");
            for i = 1:length(x)
                i != length(x) ? println(Positions_File, "$(x[i]),\t$(y[i]),\t$(z[i]),") : println(Positions_File, "$(x[i]),\t$(y[i]),\t$(z[i])")
            end
            close(Positions_File)
            N_Image += 1
        end

        ####################################################    DECISION MAKING   ################################################################
        RN = rand(1:3);

            ####################################################    MOLECULE MOVEMENT   ################################################################
            if RN == 1 && length(x) > 1
                N_Movement += 1;
                N_Displacement += 1;
                Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted = Movement(Patch_Radius, h, L, Beta, Max_Displacement_Z, Max_Displacement, σ_p, λ_p, σ_w, λ_w, Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z)
            end

            ####################################################    MOLECULE INSERTION   ################################################################
            if RN == 2
                N_Insertion += 1;
                Pc_Random, Pc_Random_Sum, Pc_Random_N, x_Insertion, y_Insertion, z_Insertion = Random_Excluded_Volume(2σ_p, h, L, σ_w, Pc_Random, Pc_Random_Sum, Pc_Random_N, x, y, z)
                #Pc_Grid, Pc_Grid_Sum, Pc_Grid_N = Grid_Excluded_Volume(2σ_p, h, L, Pc_Grid, Pc_Grid_Sum, Pc_Grid_N, x, y, z)
                #Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N = Analytic_Excluded_Volume(2σ_p, L, V, Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N, x, y, z)
                if length(x_Insertion) > 0
                    Energy, N_Insertion_Accepted, N_Insertion_Rejected = Insertion_Mezei(Patch_Radius, h, Beta, ChemPot, L, V, R_Cut, σ_p, λ_p, σ_w, λ_w, Energy, N_Insertion_Accepted, N_Insertion_Rejected, x, y, z, x_Insertion, y_Insertion, z_Insertion, Pc_Random)
                else
                    Energy, N_Insertion_Accepted, N_Insertion_Rejected = Insertion_Metropolis(Patch_Radius, h, L, V, ChemPot, Beta, R_Cut, σ_p, λ_p, σ_w, λ_w, x, y, z, Energy, N_Insertion_Accepted, N_Insertion_Rejected)
                end
            end

            ####################################################    MOLECULE REMOVAL   ################################################################
            if RN == 3 && length(x) > 1
                N_Removal += 1;
                if length(Pc_Random) == 1
                    Pc_Interpolation = Pc_Random(collect(keys(Pc_Random))[1])
                else
                    if haskey(Pc_Random, length(x) - 1)
                        Pc_Interpolation = Pc_Random[length(x) - 1]
                    else
                        Pc_Interpolation = Interpolation(Pc_Random, length(x))
                    end
                end
                if rand() > (1 - Pc_Interpolation)^1000
                    Energy, N_Removal_Accepted, N_Removal_Rejected = Removal_Mezei(Patch_Radius, Pc_Interpolation, h, L, V, Beta, ChemPot, R_Cut, σ_p, λ_p, σ_w, λ_w, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z);
                else
                    Energy, N_Removal_Accepted, N_Removal_Rejected = Removal_Metropolis(Patch_Radius, h, L, V, Beta, ChemPot, R_Cut, σ_p, λ_p, σ_w, λ_w, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z)
                end
            end

        ####################################################    MEASUREMENT SECTION   ################################################################
        if i % MC_Measurement == 0
            if i > MC_Relaxation_Steps
                N_Measurements += 1;
                Energy_Array[N_Measurements] = Energy / length(x);
                Density_Array[N_Measurements] = length(x) / V;
                if N_Measurements > 1
                    Std_Energy[N_Measurements] = std(Energy_Array[1:N_Measurements]);
                    Std_Density[N_Measurements] = std(Density_Array[1:N_Measurements]);
                end
                Mean_Energy_Array[N_Measurements] = mean(Energy_Array[1:N_Measurements]);
                Mean_Density_Array[N_Measurements] = mean(Density_Array[1:N_Measurements]);
                g_x[N_Measurements, :] = Distribution(N_Bins, L, x);
                g_y[N_Measurements, :] = Distribution(N_Bins, L, y);
                g_z[N_Measurements, :] = Distribution(N_Bins, h + 2σ_w, z);
                PotentialFunction += Potential(Patch_Radius, N_Bins, L, h, σ_w, λ_w, x, y, z)
            end
            if N_Displacement_Accepted / N_Displacement > 0.55
                Max_Displacement *= 1.05
                Max_Displacement_Z *= 1.05
            else
                Max_Displacement *= 0.95
                Max_Displacement_Z *= 0.95
            end
            Max_Displacement < 0.05 ? Max_Displacement = 0.05 : nothing
            Max_Displacement > L / 4. ? Max_Displacement = L / 4. : nothing
            Max_Displacement_Z < 0.05 ? Max_Displacement_Z = 0.05 : nothing
            Max_Displacement_Z > (h + 2σ_w) / 8. ? Max_Displacement_Z = (h + 2σ_w) / 8. : nothing
            N_Displacement, N_Displacement_Accepted = 0, 0;
        end
    end
    close(Acceptance_File)

    ####################################################   SIMULATION ENDED, OUTPUT BEGINS   ################################################################

        ####################################################    ENERGY AND DENSITY OUTPUT FILES   ################################################################
        Energy_File = open("$Output_Route/Energy.dat", "w");
        Density_File = open("$Output_Route/Density.dat", "w");
        Mean_Energy_File = open("$Output_Route/Average_Energy.dat", "w");
        Mean_Density_File = open("$Output_Route/Average_Density.dat", "w");
        println(Energy_File, "N\tEnergy")
        println(Density_File, "N\tDensity")
        println(Mean_Energy_File, "N\tAverageEnergy\tStandardDeviation")
        println(Mean_Density_File, "N\tAverageDensity\tStandardDeviation")
        for i = 1:N_Measurements
            println(Energy_File, "$i\t$(round(Energy_Array[i], digits = 8))")
            println(Density_File, "$i\t$(round(Density_Array[N_Measurements], digits = 8))")
            println(Mean_Energy_File, "$i\t$(round(Mean_Energy_Array[i], digits = 8))\t$(round(Std_Energy[i], digits = 8))")
            println(Mean_Density_File, "$i\t$(round(Mean_Density_Array[i], digits = 8))\t$(round(Std_Density[i], digits = 8))")
        end
        close(Energy_File)
        close(Density_File)
        close(Mean_Energy_File)
        close(Mean_Density_File)
    
        ####################################################    DISTRIBUTION PROFILES AND SLAB POTENTIAL  ################################################################
        Delta_xy, Delta_z = L / N_Bins, (h + 2σ_w) / N_Bins;
        r_xy, r_z = zeros(Float64, N_Bins), zeros(Float64, N_Bins);

        g_x *= (N_Bins / V);
        g_x_Mean, g_x_Std = mean(g_x, dims = 1)', std(g_x, dims = 1)';
        g_x_Normalized = g_x ./ Bulk_Density;
        g_x_Normalized_Mean, g_x_Normalized_Std = mean(g_x_Normalized, dims = 1)', std(g_x_Normalized, dims = 1)';

        g_y *= (N_Bins / V);
        g_y_Mean, g_y_Std = mean(g_y, dims = 1)', std(g_y, dims = 1)';
        g_y_Normalized = g_y ./ Bulk_Density;
        g_y_Normalized_Mean, g_y_Normalized_Std = mean(g_y_Normalized, dims = 1)', std(g_y_Normalized, dims = 1)';

        g_z *= (N_Bins / V);
        g_z_Mean, g_z_Std = mean(g_z, dims = 1)', std(g_z, dims = 1)';
        g_z_Normalized = g_z ./ Bulk_Density; 
        g_z_Normalized_Mean, g_z_Normalized_Std = mean(g_z_Normalized, dims = 1)', std(g_z_Normalized, dims = 1)';
        g_z_Max = findmax(g_z_Mean)[1] + g_z_Std[findmax(g_z_Mean)[2][1]];
        g_z_Normalized_Max = findmax(g_z_Normalized_Mean)[1] + g_z_Normalized_Std[findmax(g_z_Normalized_Mean)[2][1]];

        PotentialFunction /= N_Measurements;

        Density_Distribution_File = open("$Output_Route/Density_Distribution.dat", "w");
        Potential_File = open("$Output_Route/Potential_Function.dat", "w")
        println(Density_Distribution_File, "r_xy\tDensity_x\tStd_x\tNormalized_Density_x\tStd_Normalized_x\tDensity_y\tStd_y\tNormalized_Density_y\tStd_Normalized_y\tr_z\tDensity_z\tStd_z\tNormalized_Density_z\tStd_Normalized_z")
        println(Potential_File, "z\tU_wall(r)")
        @inbounds for i = 1:N_Bins
            r_xy[i] = round( - L / 2 + (i - 0.5) * Delta_xy, digits = 6);
            r_z[i] = round( - (h + 2σ_w) / 2 + (i - 0.5) * Delta_z, digits = 6);
            println(Density_Distribution_File, "$(r_xy[i])\t$(round(g_x_Mean[i], digits = 6))\t$(round(g_x_Std[i], digits = 6))\t$(round(g_x_Normalized_Mean[i], digits = 6))\t$(round(g_x_Normalized_Std[i], digits = 6))\t$(round(g_y_Mean[i], digits = 6))\t$(round(g_y_Std[i], digits = 6))\t$(round(g_y_Normalized_Mean[i], digits = 6))\t$(round(g_y_Normalized_Std[i], digits = 6))\t$(r_z[i])\t$(round(g_z_Mean[i], digits = 6))\t$(round(g_z_Std[i], digits = 6))\t$(round(g_z_Normalized_Mean[i], digits = 6))\t$(round(g_z_Normalized_Std[i], digits = 6))")
            println(Potential_File, "$(r_z[i])\t$(round(PotentialFunction[i], digits = 6))")
        end
        close(Density_Distribution_File)
        close(Potential_File)

        Distribution_X_Plot = Distributions_Plots("Density_Distribution_X", r_xy, g_x_Mean, g_x_Std, g_z_Max, h, Output_Route, true)
        Distribution_Y_Plot = Distributions_Plots("Density_Distribution_Y", r_xy, g_y_Mean, g_y_Std, g_z_Max, h, Output_Route, true)
        Distribution_Z_Plot = Distributions_Plots("Density_Distribution_Z", r_z, g_z_Mean, g_z_Std, g_z_Max, h, Output_Route, true)
        Distribution_X_Plot, Distribution_Y_Plot, Distribution_Z_Plot = Distribution_Plot_Unified("Density_Distribution", "Density", Distribution_X_Plot, Distribution_Y_Plot, Distribution_Z_Plot, Output_Route)

        Normalized_Distribution_X_Plot = Distributions_Plots("Normalized_Density_Distribution_X", r_xy, g_x_Mean, g_x_Std, g_z_Max, h, Output_Route, false)
        Normalized_Distribution_Y_Plot = Distributions_Plots("Normalized_Density_Distribution_Y", r_xy, g_y_Mean, g_y_Std, g_z_Max, h, Output_Route, false)
        Normalized_Distribution_Z_Plot = Distributions_Plots("Normalized_Density_Distribution_Z", r_z, g_z_Mean, g_z_Std, g_z_Max, h, Output_Route, false)
        Normalized_Distribution_X_Plot, Normalized_Distribution_Y_Plot, Normalized_Distribution_Z_Plot = Distribution_Plot_Unified("Normalized_Density_Distribution", "Normalized Density", Normalized_Distribution_X_Plot, Normalized_Distribution_Y_Plot, Normalized_Distribution_Z_Plot, Output_Route)

        Potential_Plot = plot(r_z, PotentialFunction, xlabel = "Z Axis", ylabel = "U_Wall(r)", title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        savefig(Potential_Plot, "$Output_Route/Potential_Plot")

        ####################################################    CAVITY PROBABILITY AND MOVEMENT/INSERTION/DELETION OF MOLECULES  ################################################################
        Pc_Random_Array = zeros(Float64, length(Pc_Random) - 1)
        #Pc_Grid_Array = zeros(Float64, length(Pc) - 1)
        #Pc_Analytic_Array = zeros(Float64, length(Pc) - 1)
        Pc_Density = zeros(Float64, length(Pc_Random_Array))
        @inbounds for i in keys(Pc_Random)
            if i != 0
                Pc_Random_Array[i] = Pc_Random[i];
                #Pc_Grid_Array[i] = Pc_Grid[i];
                #Pc_Analytic_Array[i] = Pc_Analytic[i];
            end
        end
        Pc_File = open("$Output_Route/Cavity_Probability.dat", "w");
        println(Pc_File, "Density\tPc_Random\tPc_Grid\tPc_Analytic")
        @inbounds for i = 1:length(Pc_Random_Array)
            Pc_Density[i] = i / V;
            println(Pc_File, "$(Pc_Density[i])\t$(round(Pc_Random_Array[i], digits = 6))")
            #println(Pc_File, "$i\t$(round(Pc_Array[i], digits = 6))\t$(round(Pc_Grid_Array[i], digits = 6))\t$(round(Pc_Analytic_Array[i], digits = 6))")
        end
        close(Pc_File)
        Cavity_Probability_Plot = plot(Pc_Density, Pc_Random_Array, xlabel = "Density", ylabel = "Cavity Probability", ylim = (0, 1), title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        savefig(Cavity_Probability_Plot, "$Output_Route/Cavity_Probability_Plot")

        Acceptance = CSV.read("$Output_Route/Acceptance.dat", delim = "\t")
        Movement_Plot = plot(Acceptance.N, Acceptance.Movement, xlabel = "Measurement", ylabel = "Movement Acceptance Percentage", ylim = (0, 100), title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        Insertion_Plot = plot(Acceptance.N, Acceptance.Insertion, xlabel = "Measurement", ylabel = "Insertion Acceptance Percentage", ylim = (0, 100), title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        Removal_Plot = plot(Acceptance.N, Acceptance.Removal, xlabel = "Measurement", ylabel = "Removal Acceptance Percentage", ylim = (0, 100), title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        savefig(Movement_Plot, "$Output_Route/Movement_Acceptance_Plot")
        savefig(Insertion_Plot, "$Output_Route/Insertion_Acceptance_Plot")
        savefig(Removal_Plot, "$Output_Route/Removal_Acceptance_Plot")

        ####################################################    ENERGY, AVERAGE ENERGY AND ENERGY HISTOGRAM  ################################################################
        println("< E / N > = $(round(Mean_Energy_Array[N_Measurements], digits = 8)) ± $(round(Std_Energy[N_Measurements], digits = 8))")

        Energy_Plot = plot(Energy_Array, xlabel = "Measurements", ylabel = "Energy [Unitless]", title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        hline!(Energy_Plot, [Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

        Energy_Histogram =  histogram(Energy_Array, xlabel = "Energy [Unitless]", ylabel = "Normalized Frequency", legend = false, framestyle = :box, bins = 20, normalize = true, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        vline!([Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
        plot!(Normal(Mean_Energy_Array[N_Measurements], Std_Energy[N_Measurements]), width = 3, linecolor = :black)
    
        Average_Energy_Plot = plot(Mean_Energy_Array, ribbon = Std_Energy, xlabel = "Measurements", ylabel = "Average Energy [Unitless]", fillalpha = 0.2, legend = false, framestyle = :box, width = 3, guidefontsize = 20, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        hline!([Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
        hline!([Mean_Energy_Array[N_Measurements] + Std_Energy[N_Measurements]], color = :black, width = 2, linestyle = :dash, linealpha = 0.5)
        hline!([Mean_Energy_Array[N_Measurements] - Std_Energy[N_Measurements]], color = :black, width = 2, linestyle = :dash, linealpha = 0.5)

        Energy_Plots = plot(Energy_Plot, Energy_Histogram, Average_Energy_Plot, layout = (@layout [a{0.3h} ; b c]))
        savefig(Energy_Plots, "$Output_Route/Energy_Plots")

        xlabel!(Energy_Plot, ""), ylabel!(Energy_Plot, "")
        xlabel!(Energy_Histogram, ""), ylabel!(Energy_Histogram, "")
        xlabel!(Average_Energy_Plot, ""), ylabel!(Average_Energy_Plot, "")

        ####################################################    DENSITY, AVERAGE DENSITY AND DENSITY HISTOGRAM  ################################################################
        println("< Density > = $(round(Mean_Density_Array[N_Measurements], digits = 8)) ± $(round(Std_Density[N_Measurements], digits = 8))")
        println("< N > = $(round(V*Mean_Density_Array[N_Measurements], digits = 8)) ± $(round(V*Std_Density[N_Measurements], digits = 8))")

        Density_Plot = plot(Density_Array, xlabel = "Measurements", ylabel = "Density [Unitless]", title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        hline!(Density_Plot, [Mean_Density_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)

        Density_Histogram =  histogram(Density_Array, xlabel = "Density [Unitless]", ylabel = "Normalized Frequency", legend = false, framestyle = :box, bins = 20, normalize = true, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        vline!([Mean_Density_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
        plot!(Normal(Mean_Density_Array[N_Measurements], Std_Density[N_Measurements]), width = 3, linecolor = :black)
    
        Average_Density_Plot = plot(Mean_Density_Array, ribbon = Std_Density, xlabel = "Measurements", ylabel = "Average Density [Unitless]", fillalpha = 0.2, legend = false, framestyle = :box, width = 3, guidefontsize = 20, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
        hline!([Mean_Density_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
        hline!([Mean_Density_Array[N_Measurements] + Std_Density[N_Measurements]], color = :black, width = 2, linestyle = :dash, linealpha = 0.5)
        hline!([Mean_Density_Array[N_Measurements] - Std_Density[N_Measurements]], color = :black, width = 2, linestyle = :dash, linealpha = 0.5)

        Density_Plots = plot(Density_Plot, Density_Histogram, Average_Density_Plot, layout = (@layout [a{0.3h} ; b c]))
        savefig(Density_Plots, "$Output_Route/Density_Plots")

        xlabel!(Density_Plot, ""), ylabel!(Density_Plot, "")
        xlabel!(Density_Histogram, ""), ylabel!(Density_Histogram, "")
        xlabel!(Average_Density_Plot, ""), ylabel!(Average_Density_Plot, "")
    
    ####################################################    SUMMARY OF THE SIMULATION  ################################################################
    Summary_File = open("$Output_Route/Summary.dat", "w")
    println(Summary_File, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   INPUT   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    println(Summary_File, "Chemical Potential = $ChemPot\nL = $L\tV = $V\nT = $T\n$MC_Relaxation_Steps Relaxation Steps.\t$MC_Relaxation_Measurement 'Measurements' every $MC_Measurement steps.\n$MC_Equilibrium_Steps Equilibrium Steps.\t$MC_Equilibrium_Measurement Measurements every $MC_Measurement steps.")
    println(Summary_File, "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   OUTPUT   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    println(Summary_File, "< E / N > = $(round(Mean_Energy_Array[N_Measurements], digits = 8)) ± $(round(Std_Energy[N_Measurements], digits = 8))\n< Density > = $(round(Mean_Density_Array[N_Measurements], digits = 8)) ± $(round(Std_Density[N_Measurements], digits = 8))\n< N > = $(round(V*Mean_Density_Array[N_Measurements], digits = 8)) ± $(round(V*Std_Density[N_Measurements], digits = 8))")
    close(Summary_File)

    ####################################################    POVRAY .pov AND .ini FILES GENERATION  ################################################################
    Povray_ini(h, L, ChemPot, T, λ_w, N_Image - 1, Patch_Percentage)
    Povray_Pov(h, L, ChemPot, T, σ_w, λ_w, Patch_Radius, Patch_Percentage)

    return Mean_Energy_Array[N_Measurements], Std_Energy[N_Measurements], Mean_Density_Array[N_Measurements], Std_Density[N_Measurements], Energy_Plot, Energy_Histogram, Average_Energy_Plot, Density_Plot, Density_Histogram, Average_Density_Plot, Distribution_X_Plot, Distribution_Y_Plot, Distribution_Z_Plot, Normalized_Distribution_X_Plot, Normalized_Distribution_Y_Plot, Normalized_Distribution_Z_Plot
end

function Movement(Patch_Radius::Float64, h::Float64, L::Float64, Beta::Float64, z_Displacement::Float64, Displacement::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, Energy::Float64, N_Movement_Accepted::Int64, N_Movement_Rejected::Int64, N_Displacement_Accepted::Int64, R_Cut::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    j = rand(1:length(x))
    Energy_Old = Energy_Calculation(Patch_Radius, h, L, R_Cut, σ_p, λ_p, σ_w, λ_w, x[j], y[j], z[j], x, y, z)
    x_Old, y_Old, z_Old = x[j], y[j], z[j];
    x[j] += Displacement * (rand() - 0.5);
    x[j] = PeriodicBoundaryConditions(L, x[j]);
    y[j] += Displacement * (rand() - 0.5);
    y[j] = PeriodicBoundaryConditions(L, y[j]);
    z[j] += z_Displacement * (rand() - 0.5);
    Energy_New = Energy_Calculation(Patch_Radius, h, L, R_Cut, σ_p, λ_p, σ_w, λ_w, x[j], y[j], z[j], x, y, z);
    if rand() < exp(-Beta * (Energy_New - Energy_Old) )
        N_Movement_Accepted += 1;
        N_Displacement_Accepted += 1;
        Energy += (Energy_New - Energy_Old);
    else
        N_Movement_Rejected += 1;
        x[j], y[j], z[j] = x_Old, y_Old, z_Old;
    end
    return Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted
end

function Insertion_Metropolis(Patch_Radius::Float64, h::Float64, L::Float64, V::Float64, ChemPot::Float64, Beta::Float64, R_Cut::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, Energy::Float64, N_Insertion_Accepted::Int64, N_Insertion_Rejected::Int64)
    x_Insertion, y_Insertion, z_Insertion = L * (rand() - 0.5), L * (rand() - 0.5), (h + 2σ_w) * (rand() - 0.5)
    Energy_Insertion = Energy_Calculation(Patch_Radius, h, L, R_Cut, σ_p, λ_p, σ_w, λ_w, x_Insertion, y_Insertion, z_Insertion, x, y, z)
    if rand() < exp( Beta * (ChemPot - Energy_Insertion) + log(V / (length(x) + 1)) )
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion), append!(y, y_Insertion), append!(z, z_Insertion)
        Energy += Energy_Insertion;
    else
        N_Insertion_Rejected += 1;
    end
    return Energy, N_Insertion_Accepted, N_Insertion_Rejected
end

function Insertion_Mezei(Patch_Radius::Float64, h::Float64, Beta::Float64, ChemPot::Float64, L::Float64, V::Float64, R_Cut::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, Energy::Float64, N_Insertion_Accepted::Int64, N_Insertion_Rejected::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, x_Insertion::Array{Float64, 1}, y_Insertion::Array{Float64, 1}, z_Insertion::Array{Float64, 1}, Pc::Dict{Int64, Float64})
    j = rand(1:length(x_Insertion))
    Energy_Insertion = Energy_Calculation(Patch_Radius, h, L, R_Cut, σ_p, λ_p, σ_w, λ_w, x_Insertion[j], y_Insertion[j], z_Insertion[j], x, y, z);
    if rand() < (V * Pc[length(x)] / (length(x) + 1) ) * exp(Beta * (ChemPot - Energy_Insertion))
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion[j]), append!(y, y_Insertion[j]), append!(z, z_Insertion[j])
        Energy += Energy_Insertion;
    else
        N_Insertion_Rejected += 1;
    end
    return Energy, N_Insertion_Accepted, N_Insertion_Rejected
end

function Removal_Metropolis(Patch_Radius::Float64, h::Float64, L::Float64, V::Float64, Beta::Float64, ChemPot::Float64, R_Cut::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, Energy::Float64, N_Removal_Accepted::Int64, N_Removal_Rejected::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    j = rand(1:length(x));
    Energy_Removal = Energy_Calculation(Patch_Radius, h, L, R_Cut, σ_p, λ_p, σ_w, λ_w, x[j], y[j], z[j], x, y, z)
    if rand() < exp( Beta * (Energy_Removal - ChemPot) + log(length(x) / V) )
        N_Removal_Accepted += 1;
        deleteat!(x, j), deleteat!(y, j), deleteat!(z, j)
        Energy -= Energy_Removal;
    else
        N_Removal_Rejected += 1;
    end
    return Energy, N_Removal_Accepted, N_Removal_Rejected 
end

function Removal_Mezei(Patch_Radius::Float64, Pc::Float64, h::Float64, L::Float64, V::Float64, Beta::Float64, ChemPot::Float64, R_Cut::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, Energy::Float64, N_Removal_Accepted::Int64, N_Removal_Rejected::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    j = rand(1:length(x))
    Energy_Removal = Energy_Calculation(Patch_Radius, h, L, R_Cut, σ_p, λ_p, σ_w, λ_w, x[j], y[j], z[j], x, y, z)
    if rand() < ( length(x) / (V * Pc) ) * exp(Beta * (Energy_Removal - ChemPot))
        N_Removal_Accepted += 1;
        deleteat!(x, j), deleteat!(y, j), deleteat!(z, j)
        Energy -= Energy_Removal;
    else
        N_Removal_Rejected += 1;
    end
    return Energy, N_Removal_Accepted, N_Removal_Rejected
end

function Energy_Calculation(Patch_Radius::Float64, h::Float64, L::Float64, R_Cut::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, rx::Float64, ry::Float64, rz::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    Energy = 0;
    ####################################################    ENERGY CONTRIBUTION FROM THE SLAB  ################################################################
    if abs(rz) > h / 2 - λ_w
        x_Index = ceil(rx + L / 2) / (2σ_w) + 2;
        y_Index = ceil((ry + L / 2) / (2 * √3 * σ_w) + √3);

        for i_y = y_Index - 2:y_Index + 2, i_x = x_Index - 2:x_Index + 2
            x_Position = - L / 2 - 3σ_w + (i_x - 1) * 2σ_w;
            y_Position = - L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w;
            if x_Position < rx + λ_w && x_Position > rx - λ_w
                    if y_Position < ry + λ_w && y_Position > ry - λ_w
                        Delta_x = rx - x_Position;
                        Delta_y = ry - y_Position;
                        Delta_z = abs(rz) - (h / 2 + σ_w)
                        r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
                        r_center = sqrt( (x_Position)^2 + (y_Position)^2)
                        r_center < Patch_Radius ? Energy += u_SquareWell(r2, σ_w, λ_w, 1.0) : Energy += u_SquareWell(r2, σ_w, λ_w, 0.0)
                        Energy == Inf ? (return Energy) : nothing
                    end
            end

            x_Position = - L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w;
            y_Position = - L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w;
            if x_Position < rx + λ_w && x_Position > rx - λ_w
                    if y_Position < ry + λ_w && y_Position > ry - λ_w
                        Delta_x = rx - x_Position;
                        Delta_y = ry - y_Position;
                        Delta_z = abs(rz) - (h / 2 + σ_w)
                        r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
                        r_center = sqrt( (x_Position)^2 + (y_Position)^2 );
                        r_center < Patch_Radius ? Energy += u_SquareWell(r2, σ_w, λ_w, 1.0) : Energy += u_SquareWell(r2, σ_w, λ_w, 0.0);
                        Energy == Inf ? (return Energy) : nothing
                    end
            end
        end
    end

    ####################################################    ENERGY CONTRIBUTION BETWEEN PARTICLES  ################################################################
    @inbounds for i = 1:length(x)
        Delta_x = rx - x[i];
        Delta_x = PeriodicBoundaryConditions(L, Delta_x);
        Delta_y = ry - y[i];
        Delta_y = PeriodicBoundaryConditions(L, Delta_y);
        Delta_z = rz - z[i];
        r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
        if r2 != 0.
            if r2 < R_Cut^2
                Energy += u_SquareWell(r2, σ_p, λ_p);
            end
        end
        Energy == Inf ? (return Energy) : nothing
    end

    return Energy
end

function u_SquareWell(r2::Float64, σ::Float64, λ::Float64, e::Float64 = 1.)
    if r2 <= (2σ)^2
        return Inf
    elseif r2 <= λ^2
        return -e
    else
        return 0
    end
end

function PeriodicBoundaryConditions(L::Float64, x::Float64)
    return x - L * round(x / L)
end

function Random_Excluded_Volume(Overlap::Float64, h::Float64, L::Float64, σ_w::Float64, Pc::Dict{Int64, Float64}, Pc_Sum::Dict{Int64, Float64}, Pc_N::Dict{Int64, Int64}, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    N_Random = 1000;
    N_in = 0;
    x_Insertion, y_Insertion, z_Insertion = Float64[], Float64[], Float64[];
    @inbounds for i = 1:N_Random
        Control = false;
        z_V = (h + 2σ_w) * (rand() - 0.5);
        x_V = L * (rand() - 0.5);
        y_V = L * (rand() - 0.5);
        for j = 1:length(x)
            Delta_x = x_V - x[j];
            Delta_x = PeriodicBoundaryConditions(L, Delta_x);
            Delta_y = y_V - y[j];
            Delta_y = PeriodicBoundaryConditions(L, Delta_y);
            Delta_z = z_V - z[j];
            r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
            if r2 < Overlap^2
                Control = true;
                break
            end
        end
        if Control == false
            for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
                Delta_x = x_V - (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w);
                Delta_y = y_V - (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w);
                Delta_z = abs(z_V) - (h / 2 + σ_w)
                r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
                if r2 < Overlap^2
                    Control = true;
                    break
                end
                Delta_x = x_V - (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w);
                Delta_y = y_V - (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w);
                Delta_z = abs(z_V) - (h / 2 + σ_w)
                r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
                if r2 < Overlap^2
                    Control = true;
                    break
                end
            end
        end
        if Control == false
            N_in += 1;
            append!(x_Insertion, x_V)
            append!(y_Insertion, y_V)
            append!(z_Insertion, z_V)
        end
    end
    Volume_Ratio = N_in /N_Random;
    if !haskey(Pc, length(x))
        Pc_N[length(x)] = 1
        Pc_Sum[length(x)] = Volume_Ratio;
        Pc[length(x)] = Volume_Ratio;
    else
        Pc_N[length(x)] += 1;
        Pc_Sum[length(x)] += Volume_Ratio;
        Pc[length(x)] = Pc_Sum[length(x)] / Pc_N[length(x)]
    end
    return Pc, Pc_Sum, Pc_N, x_Insertion, y_Insertion, z_Insertion
end

function Grid_Excluded_Volume(Overlap::Float64, h::Float64, L::Float64, Pc_Grid::Dict{Int64, Float64}, Pc_Grid_Sum::Dict{Int64, Float64}, Pc_Grid_N::Dict{Int64, Int64}, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    Grid = 10;
    Delta = L / (Grid + 1);
    Δz = h / (Grid + 1)
    EVMPS = zeros(Grid, Grid, Grid);
    @inbounds for i_x = 1:Grid
        @inbounds for i_y = 1:Grid
            @inbounds for i_z = 1:Grid
                x_Grid = i_x * Delta - L / 2.;
                y_Grid = i_y * Delta - L / 2.;
                z_Grid = i_z * Δz - h / 2.;
                @inbounds for i = 1:length(x)
                    Delta_x = x_Grid - x[i];
                    Delta_y = y_Grid - y[i];
                    Delta_z = z_Grid - z[i];
                    r2 = Delta_x^2 + Delta_y^2 + Delta_z^2
                    if r2 < Overlap^2
                        EVMPS[i_x, i_y, i_z] = 1;
                        break
                    end
                end
            end
        end
    end
    Volume_Ratio = 1 - sum(EVMPS) / Grid ^ 3;
    if !haskey(Pc_Grid, length(x))
        Pc_Grid_N[length(x)] = 1;
        Pc_Grid_Sum[length(x)] = Volume_Ratio;
        Pc_Grid[length(x)] = Volume_Ratio;
    else
        Pc_Grid_N[length(x)] += 1;
        Pc_Grid_Sum[length(x)] += Volume_Ratio;
        Pc_Grid[length(x)] = Pc_Grid_Sum[length(x)] / Pc_Grid_N[length(x)];
    end
    return Pc_Grid, Pc_Grid_Sum, Pc_Grid_N
end

function Analytic_Excluded_Volume(Overlap::Float64, L::Float64, V::Float64, Pc_Analytic::Dict{Int64, Float64}, Pc_Analytic_Sum::Dict{Int64, Float64}, Pc_Analytic_N::Dict{Int64, Int64}, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    V_Excluded = length(x) * (4. / 3.) * pi * Overlap^3;
    V_Excluded_Correction = 0;
    @inbounds for i = 1:length(x) - 1
        @inbounds for j = i + 1:length(x)
            Delta_x = x[i] - x[j];
            Delta_y = y[i] - y[j];
            Delta_z = z[i] - z[j];
            r2 = Delta_x^2. + Delta_y^2. + Delta_z^2.;
            if r2 < (2Overlap)^2.
                d = sqrt(r2);
                V_Excluded_Correction += (pi / 12.) * (4 * Overlap + d) * (2 * Overlap - d)^2.;
            else
                Delta_x > L - 1 ? Delta_x -= L : nothing
                Delta_y > L - 1 ? Delta_y -= L : nothing
                r2 = Delta_x^2. + Delta_y^2. + Delta_z^2.
                if r2 < (2Overlap)^2.
                    d = sqrt(r2);
                    V_Excluded_Correction += (pi / 12.) * (4 * Overlap + d) * (2 * Overlap - d)^2.;
                end
            end
        end
    end
    Volume_Ratio = 1 - (V_Excluded - V_Excluded_Correction) / V;
    Volume_Ratio > 1 || Volume_Ratio < 0 ? error("Volume Ratio ($Volume_Ratio) can't be negative or greater than one.") : nothing
    if !haskey(Pc_Analytic, length(x))
        Pc_Analytic_N[length(x)] = 1;
        Pc_Analytic_Sum[length(x)] = Volume_Ratio;
        Pc_Analytic[length(x)] = Volume_Ratio;
    else
        Pc_Analytic_N[length(x)] += 1;
        Pc_Analytic_Sum[length(x)] += Volume_Ratio;
        Pc_Analytic[length(x)] = Pc_Analytic_Sum[length(x)] / Pc_Analytic_N[length(x)];
    end
    return Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N
end

function Interpolation(Pc::Dict{Int64, Float64}, l::Float64)
    if length(Pc) == 1
        Pc_Interpolation = collect(keys(Pc))[1];
    elseif length(Pc) > 1
        lim_inf = l -1
        while true
            if haskey(Pc, lim_inf)
                break
            end
            lim_inf -= 1;
        end

        lim_sup = l + 1;
        while true
            if haskey(Pc, lim_sup)
                break
            end
            lim_sup += 1;
        end
    elseif length(Pc) == 0
        error("Pc has no values still.")
    end
    Pc_Interpolation = Pc[lim_inf] * ((1 - l - lim_inf) / (lim_sup - lim_inf)) + Pc[lim_sup] * ((l - lim_inf) / (lim_sup - lim_inf))
    return Pc_Interpolation
end

function Potential(Patch_Radius::Float64, N_Bins::Int64, L::Float64, h::Float64, σ_w::Float64, λ_w::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    Delta = (h + 2σ_w) / N_Bins;
    PotentialFunction = zeros(Float64, N_Bins);
    N_Slab = 2 * convert(Int64, ceil(L / (2σ_w) + 3 + 1) * ceil(L / (√3 * 2σ_w) + √3))
    for i = 1:length(x)
        l = convert(Int64, ceil( ((h + 2σ_w) / 2. - z[i]) / Delta) )
        Energy = 0;
        for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
            Delta_x = x[i] - (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w);
            Delta_y = y[i] - (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w);
            Delta_z = abs(z[i]) - (h / 2 + σ_w)
            r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
            r_center = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)^2 + (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)^2)
            if r_center < Patch_Radius
                Energy += u_SquareWell(r2, σ_w, λ_w, 1.0);
            else
                Energy += u_SquareWell(r2, σ_w, λ_w, 0.0);
            end
            Delta_x = x[i] - (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w);
            Delta_y = y[i] - (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w);
            Delta_z = abs(z[i]) - (h / 2 + σ_w)
            r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
            r_center = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)^2 + ((- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w))^2 );
            if r_center < Patch_Radius
                Energy += u_SquareWell(r2, σ_w, λ_w, 1.0);
            else
                Energy += u_SquareWell(r2, σ_w, λ_w, 0.0);
            end
            Energy == Inf ? error("Particle inside infinite potential") : nothing
        end
        PotentialFunction[l] += Energy;
    end
    return PotentialFunction
end

function Distribution(N_Bins::Int64, L::Float64, x::Array{Float64, 1})
    Delta = L / N_Bins;
    g = zeros(Float64, N_Bins);
    for i = 1:length(x)
        l = convert(Int64, ceil( (L / 2. - x[i]) / Delta) )
        g[l] += 1
    end
    return g
end

function Povray_Pov(h::Float64, L::Float64, ChemPot::Float64, T::Float64, σ_w::Float64, λ_w::Float64, Patch_Radius::Float64, Patch_Percentage::Int64)
    Particle_r, Particle_g, Particle_b, Particle_t = 174/255, 214/255, 241/255, 0;
    Slab_Attractive_r, Slab_Attractive_g, Slab_Attractive_b, Slab_Attractive_t = 241/255, 196/255, 15/255, 0; 
    Slab_Neutral_r, Slab_Neutral_g, Slab_Neutral_b, Slab_Neutral_t  = 142/255, 68/255, 173/255, 0; 
    Slab_Repulsive_r, Slab_Repulsive_g, Slab_Repulsive_b, Slab_Repulsive_t = 183/255, 28/255, 28/255, 0;
    Wall_r, Wall_g, Wall_b, Wall_t = 174/255, 214/255, 241/255, 0.5;
    Output_Route = pwd() * "/Output_UPDATING/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/h_$(h)/Positions"
    ##################################################      X AXIS VIEW     ################################################################
    Pov_File = open("$Output_Route/Pore_X_Axis.pov", "w");
    println(Pov_File, "global_settings {\n\tambient_light rgb <0.2, 0.2, 0.2>\tmax_trace_level 15\n}\n")
    println(Pov_File, "background { color rgb <1, 1, 1> }\n")
    println(Pov_File, "#default { finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.7 phong 1} }\n")
    h > L ? println(Pov_File, "camera {\n\tperspective\n\tlocation <0, $(-1.5h), 0>\n\tlook_at <0, 0, 0>\n}\n") : println(Pov_File, "camera {\n\tperspective\n\tlocation <0, $(-1.5L), 0>\n\tlook_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(-5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(+5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(-5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(+5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "#macro Particle(rx, ry, rz)\n\tintersection {\n\t\t\tsphere {\n\t\t\t<rx, ry, rz>, 0.5\n\t\t\tpigment {rgbt <$Particle_r, $Particle_g, $Particle_b, $Particle_t> }\n\t\t}\n\t\tbox {\n\t\t\t<-L/2, -L/2, h/2>,\t<L/2, L/2, -h/2>\n\t\t\tpigment {rgbt <$Particle_r, $Particle_g, $Particle_b, $Particle_t> }\n\t\t}\n\tno_shadow}\n#end\n")
    println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t#if (sqrt(pow(rx, 2) + pow(ry, 2)) <= $Patch_Radius)\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\nsphere {\n\t\t<rx, ry, rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n#else\nsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Neutral_r, $Slab_Neutral_g, $Slab_Neutral_b, $Slab_Neutral_t> }\n\tno_shadow}\nsphere {\n\t\t<rx, ry, rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Neutral_r, $Slab_Neutral_g, $Slab_Neutral_b, $Slab_Neutral_t> }\n\tno_shadow}\n\t#end\n#end")
    println(Pov_File, "#macro Wall(L, h)\n\tunion{\n\t\ttriangle {\n\t\t\t<-L / 2, -L / 2, h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t\ttriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t}
        union {\n\t\ttriangle {\n\t\t\t<L / 2, -L / 2, h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n\n\t
        union {\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <-L / 2, L / 2, -h /2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\ntriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <L / 2, L / 2, h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n#end")
    println(Pov_File, "#declare L = $L;\n#declare h = $h;")
    println(Pov_File, """#fopen File_Positions concat("$Output_Route/Pos_", str(clock, 1, 0), ".xyz") read""")
    println(Pov_File, "\t#while (defined( File_Positions ))\n\t\t#read (File_Positions, rx, ry, rz)\n\t\tParticle(rx, ry, -rz)\n\t\t#declare PBC = false;\n\t\t#if (rx > (L - 1) / 2)\n\t\t\t#declare rx = rx - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (rx < -(L - 1) / 2)\n\t\t\t#declare rx = rx + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t
        #if (ry > (L - 1) / 2)\n\t\t\t#declare ry = ry - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (ry < -(L - 1) / 2)\n\t\t\t#declare ry = ry + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (PBC)\n\t\t\tParticle(rx, ry, -rz)\n\t\t#end\n\t#end\n#fclose File_Positions")
    println(Pov_File, """#fopen File_Slab "$Output_Route/Slab.xyz" read""")
    println(Pov_File, "#while (defined (File_Slab)) \n\t#read (File_Slab, rx, ry, rz)\n\tSlab(rx, ry, rz, $σ_w, $λ_w)\n#end\n#fclose File_Slab")
    println(Pov_File, "Wall(L, h)")
    close(Pov_File)
    ##################################################      Z AXIS VIEW     ################################################################
    Pov_File = open("$Output_Route/Pore_Z_Axis.pov", "w");
    println(Pov_File, "global_settings {\n\tambient_light rgb <0.2, 0.2, 0.2>\tmax_trace_level 15\n}\n")
    println(Pov_File, "background { color rgb <1, 1, 1> }\n")
    println(Pov_File, "#default { finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.7 phong 1} }\n")
    println(Pov_File, "camera {\n\tperspective\n\tlocation <0, 0, $(1.3L)>\n\tlook_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(-5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(+5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(-5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(+5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "#macro Particle(rx, ry, rz)\n\tintersection {\n\t\t\tsphere {\n\t\t\t<rx, ry, rz>, 0.5\n\t\t\tpigment {rgbt <$Particle_r, $Particle_g, $Particle_b, $Particle_t> }\n\t\t}\n\t\tbox {\n\t\t\t<-L/2, -L/2, h/2>,\t<L/2, L/2, -h/2>\n\t\t\tpigment {rgbt <$Particle_r, $Particle_g, $Particle_b, $Particle_t> }\n\t\t}\n\tno_shadow}\n#end\n")
    println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t#if (sqrt(pow(rx, 2) + pow(ry, 2)) <= $Patch_Radius)\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n#else\nsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Neutral_r, $Slab_Neutral_g, $Slab_Neutral_b, $Slab_Neutral_t> }\n\tno_shadow}\n\t#end\n#end")
    println(Pov_File, "#macro Wall(L, h)\n\tunion{\n\t\ttriangle {\n\t\t\t<-L / 2, -L / 2, h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t\ttriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t}
        union {\n\t\ttriangle {\n\t\t\t<L / 2, -L / 2, h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n\n\t
        union {\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <-L / 2, L / 2, -h /2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\ntriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <L / 2, L / 2, h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n
        union {\n\t\ttriangle {\n\t\t\t<L / 2, -L / 2, -h /2>, <L / 2, -L / 2, h / 2>, <-L / 2, -L / 2, -h /2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\ntriangle {\n\t\t\t<-L / 2, -L / 2, -h /2>, <-L / 2, -L / 2, h / 2>, <L / 2, -L / 2, h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n#end")
    println(Pov_File, "#declare L = $L;\n#declare h = $h;")
    println(Pov_File, """#fopen File_Positions concat("$Output_Route/Pos_", str(clock, 1, 0), ".xyz") read""")
    println(Pov_File, "\t#while (defined( File_Positions ))\n\t\t#read (File_Positions, rx, ry, rz)\n\t\tParticle(rx, ry, -rz)\n\t\t#declare PBC = false;\n\t\t#if (rx > (L - 1) / 2)\n\t\t\t#declare rx = rx - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (rx < -(L - 1) / 2)\n\t\t\t#declare rx = rx + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t
        #if (ry > (L - 1) / 2)\n\t\t\t#declare ry = ry - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (ry < -(L - 1) / 2)\n\t\t\t#declare ry = ry + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (PBC)\n\t\t\tParticle(rx, ry, -rz)\n\t\t#end\n\t#end\n#fclose File_Positions")
    println(Pov_File, """#fopen File_Slab "$Output_Route/Slab.xyz" read""")
    println(Pov_File, "#while (defined (File_Slab)) \n\t#read (File_Slab, rx, ry, rz)\n\tSlab(rx, ry, rz, $σ_w, $λ_w)\n#end\n#fclose File_Slab")
    println(Pov_File, "Wall(L, h)")
    close(Pov_File)
    ##################################################      WALL DECORATION     ################################################################
    Pov_File = open("$Output_Route/Pore_Slab.pov", "w");
    println(Pov_File, "global_settings {\n\tambient_light rgb <0.2, 0.2, 0.2>\tmax_trace_level 15\n}\n")
    println(Pov_File, "background { color rgb <1, 1, 1> }\n")
    println(Pov_File, "#default { finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.7 phong 1} }\n")
    println(Pov_File, "camera {\n\tperspective\n\tlocation <0, 0, $(1.3L)>\n\tlook_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(-5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(+5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(-5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(+5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t#if (sqrt(pow(rx, 2) + pow(ry, 2)) <= $Patch_Radius)\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n#else\nsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Neutral_r, $Slab_Neutral_g, $Slab_Neutral_b, $Slab_Neutral_t> }\n\tno_shadow}\n\t#end\n#end")
    println(Pov_File, """#fopen File_Slab "$Output_Route/Slab.xyz" read""")
    println(Pov_File, "#while (defined (File_Slab)) \n\t#read (File_Slab, rx, ry, rz)\n\tSlab(rx, ry, 0, $σ_w, $λ_w)\n#end\n#fclose File_Slab")
    close(Pov_File)
end

function Povray_ini(h::Float64, L::Float64, ChemPot::Float64, T::Float64, λ_w::Float64, Frames::Int64, Patch_Percentage::Int64)
    Output_Route = pwd() * "/Output_UPDATING/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/h_$(h)/Positions"
    mkpath("$Output_Route")
    Ini_File = open("$Output_Route/Pore_X_Axis.ini", "w");
    println(Ini_File, "Input_File_Name = $Output_Route/Pore_X_Axis.pov")
    println(Ini_File, "Output_File_Name = $Output_Route/")
    println(Ini_File, "+W800 +H800\n")
    println(Ini_File, "Initial_Frame = 1")
    println(Ini_File, "Final_Frame = $Frames")
    println(Ini_File, "Initial_Clock = 1")
    println(Ini_File, "Final_Clock = $Frames\n")
    println(Ini_File, "Cyclic_Animation = off")
    close(Ini_File)

    Ini_File = open("$Output_Route/Pore_Z_Axis.ini", "w");
    println(Ini_File, "Input_File_Name = $Output_Route/Pore_Z_Axis.pov")
    println(Ini_File, "Output_File_Name = $Output_Route/")
    println(Ini_File, "+W800 +H800\n")
    println(Ini_File, "Initial_Frame = 1")
    println(Ini_File, "Final_Frame = $Frames")
    println(Ini_File, "Initial_Clock = 1")
    println(Ini_File, "Final_Clock = $Frames\n")
    println(Ini_File, "Cyclic_Animation = off")
    close(Ini_File)
end

function Pore_Plots(Plot_Name, X_Axis_Label, Y_Axis_Label, p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, Output_Route, T, ChemPot, Bulk_Density, Patch_Percentage)
    y = ones(3);
    Title = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, y[2], Plots.text("T = $T, \\mu = $(round(ChemPot, digits = 3)), \\rho_{Bulk} = $Bulk_Density", :black, 50)), axis=false, grid = false, leg=false)
    Y_Axis = plot(guide_position = :right, left_margin = -30mm, right_margin = 20mm, ylabel = "$Y_Axis_Label", guidefontsize = 30, grid = false, axis = false)
    X_Axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("$X_Axis_Label", :black, 30)), axis=false, grid = false, leg=false)
    Pore_Plot = plot(Y_Axis, Title, p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, X_Axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Pore_Plot, "$(Output_Route)$(Plot_Name)_T_$(T)_ChemPot_$(round(ChemPot, digits = 2))_Patch_$Patch_Percentage%")
end

function Distributions_Plots(Plot_Name::String, X_Array, Y_Array, Yerror_Array, Y_Max::Float64, h::Float64, Output_Route::String, Normalized::Bool)
    Distribution_Plot = plot(X_Array, Y_Array, ribbon = Yerror_Array, ylim = (0, Y_Max), legend = false, fillalpha = 0.2, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, title = "Slit Separation: $h", titlefontsize = 25, width = 3)
    Normalized ? hline!(Distribution_Plot, [1.], width = 3, linecolor = :black, linestyle = :dash) : nothing
    return Distribution_Plot
end

function Distribution_Plot_Unified(Plot_Name::String, Y_Axis_Label::String, X_Plot, Y_Plot, Z_Plot, Output_Route::String)
    xlabel!(X_Plot, "X Axis")
    xlabel!(Y_Plot, "Y Axis")
    xlabel!(Z_Plot, "Z Axis")
    ylabel!(X_Plot, "$Y_Axis_Label")
    title!(X_Plot, "")
    title!(Z_Plot, "")
    Plots_Unified = plot(X_Plot, Y_Plot, Z_Plot, layout = (1, 3), link = :y, guidefontsize = 20, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    savefig(Plots_Unified, "$(Output_Route)/$(Plot_Name)")
    xlabel!(X_Plot, "")
    xlabel!(Y_Plot, "")
    xlabel!(Z_Plot, "")
    ylabel!(X_Plot, "")
    title!(X_Plot, "Slit Separation: $h")
    title!(Z_Plot, "Slit Separation: $h")
    return X_Plot, Y_Plot, Z_Plot
end

@time Pore_Separation()