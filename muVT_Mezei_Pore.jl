using Statistics;
using Plots; gr();
using Plots.PlotMeasures;
using Test;
using CSV;
using StatsPlots;
using Distributions;

function Mezei(ChemPot::Float64, h::Float64, L::Float64, T::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, Number_Run::Int64, Total_Run::Int64, Patch_Percentage::Int64, Bulk_Density::Float64, R_Cut::Float64 = 3.)
    V = (h + 2σ_w) * L^2;
    N_Bins = 200;
    MC_Measurement = convert(Int64, ceil( 8V ));
    MC_Relaxation_Measurement = 1_000;
    MC_Equilibrium_Measurement = 10_000;
    """     CONFIGURATIONAL STEPS       """
    MC_Relaxation_Steps = MC_Measurement * MC_Relaxation_Measurement;
    MC_Equilibrium_Steps = MC_Measurement * MC_Equilibrium_Measurement;
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps;
    """     VARIABLE INITIALIZATION     """
    x, y, z = Float64[], Float64[], Float64[];
    Beta = 1. / T;
    Pc, Pc_Sum, Pc_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    #Pc_Grid, Pc_Grid_Sum, Pc_Grid_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    #Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    z_Displacement, Displacement, N_Displacement, N_Displacement_Accepted = (h + 2σ_w) / 8., L / 8., 0, 0;
    N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
    N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
    N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
    Energy, N_Measurements = 0., 0;
    Energy_Sum, Density_Sum = 0., 0.;
    Energy_Array, Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    Average_Energy_Array, Average_Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    σ_Energy, σ_Density = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ); 
    g_x, g_y, g_z = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ), N_Bins), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ), N_Bins), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ), N_Bins);
    PotentialFunction = zeros(Float64, N_Bins);
    Delta = h / N_Bins;
    ρ_Max = 1 / ( (4. / 3.) * π * σ_p^3 );
    N_Image = 1;
    Patch_Radius = round(sqrt(L^2 * Patch_Percentage / 314) , digits = 6)
    """     OUTPUT FILES        """
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/h_$(h)"
    mkpath("$Output_Route/Positions")
    Average_Energy_File = open("$Output_Route/Average_Energy.dat", "w");
    println(Average_Energy_File, "N\t< E / N >")
    Energy_File = open("$Output_Route/Energy.dat", "w");
    println(Energy_File, "N\tE / N ")
    Average_Density_File = open("$Output_Route/Average_Density.dat", "w");
    println(Average_Density_File, "N\t< Density >")
    Density_File = open("$Output_Route/Density.dat", "w");
    println(Density_File, "N\tDensity")
    Displacement_File = open("$Output_Route/DisplacementAcceptance.dat", "w");
    println(Displacement_File, "N\tAcceptance")
    Insertion_File = open("$Output_Route/InsertionAcceptance.dat", "w");
    println(Insertion_File, "N\tAcceptance")
    Removal_File = open("$Output_Route/RemovalAcceptance.dat", "w");
    println(Removal_File, "N\tAcceptance")
    """          SLAB CREATION    """
    N_Slab = 2 * convert(Int64, ceil(L / (2σ_w) + 3 + 1) * ceil(L / (√3 * 2σ_w) + √3))
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
    """     SIMULATIONS CYCLES      """
    @inbounds for i = 1:MC_Steps
        """     PRINTS PROGRESS TO SCREEN   """
        if i < MC_Relaxation_Steps && i % floor(0.01MC_Relaxation_Steps) == 0
            @test all(Array(x) .<= L / 2.) && all(Array(x) .>= -L / 2.)
            @test all(Array(y) .<= L / 2.) && all(Array(y) .>= -L / 2.)
            @test all(Array(z) .<= h / 2. + σ_w) && all(Array(z) .>= -h / 2. - σ_w)
            @test length(x) / V <= ρ_Max

            println("$(convert(Int64, floor(100i / MC_Relaxation_Steps)))% Relaxation. [$Number_Run / $Total_Run]")
            println("μ = $ChemPot\tT = $T\tρ_Bulk = $Bulk_Density")
            println("Steps: $i / $MC_Relaxation_Steps")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("N = $(length(x))")
            println("Density = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Max Z Displacement = $(round(z_Displacement, digits = 6))")
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

            println(Displacement_File, "$(convert(Int64, floor(100i / MC_Relaxation_Steps)))\t$(round(100N_Displacement_Accepted / N_Displacement, digits = 2))")
            println(Insertion_File, "$(convert(Int64, floor(100i / MC_Relaxation_Steps)))\t$(round(100N_Insertion_Accepted / N_Insertion, digits = 2))")
            println(Removal_File, "$(convert(Int64, floor(100i / MC_Relaxation_Steps)))\t$(round(100N_Removal_Accepted / N_Removal, digits = 2))")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
        end

        if i > MC_Relaxation_Steps && i % floor(0.01MC_Equilibrium_Steps) == 0
            @test all(Array(x) .<= L / 2.) && all(Array(x) .>= -L / 2.)
            @test all(Array(y) .<= L / 2.) && all(Array(y) .>= -L / 2.)
            @test all(Array(z) .<= h / 2. + σ_w) && all(Array(z) .>= -h / 2. - σ_w)
            @test length(x) / V <= ρ_Max

            println("$(convert(Int64, floor(100(i - MC_Relaxation_Steps) / MC_Equilibrium_Steps)))% Equilibrium. [$Number_Run / $Total_Run]")
            println("μ = $ChemPot\tT = $T\tρ_Bulk = $Bulk_Density")
            println("Steps: $(i - MC_Relaxation_Steps) / $MC_Equilibrium_Steps.\t[$N_Measurements Measurements]")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("N = $(length(x))")
            println("Density = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))")
            println("Max Z Displacement = $(round(z_Displacement, digits = 6))")
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
            println(Displacement_File, "$(convert(Int64, floor(100 + 100i / MC_Relaxation_Steps)))\t$(round(100N_Displacement_Accepted / N_Displacement, digits = 2))")
            println(Insertion_File, "$(convert(Int64, floor(100 + 100i / MC_Relaxation_Steps)))\t$(round(100N_Insertion_Accepted / N_Insertion, digits = 2))")
            println(Removal_File, "$(convert(Int64, floor(100 + 100i / MC_Relaxation_Steps)))\t$(round(100N_Removal_Accepted / N_Removal, digits = 2))")
            N_Movement, N_Movement_Accepted, N_Movement_Rejected = 0, 0, 0;
            N_Insertion, N_Insertion_Accepted, N_Insertion_Rejected = 0, 0, 0;
            N_Removal, N_Removal_Accepted, N_Removal_Rejected = 0, 0, 0;
        end

        if i > MC_Relaxation_Steps && i % .1MC_Equilibrium_Steps == 0
            Positions_File = open("$Output_Route/Positions/Pos_$N_Image.xyz", "w");
            for i = 1:length(x)
                i != length(x) ? println(Positions_File, "$(x[i]),\t$(y[i]),\t$(z[i]),") : println(Positions_File, "$(x[i]),\t$(y[i]),\t$(z[i])")
            end
            close(Positions_File)
            N_Image += 1
        end

        i == MC_Relaxation_Steps ? println("~~~    STARTING MEASUREMENT STEPS    ~~~") : nothing
        RN = rand(1:3);
        if RN == 1 && length(x) > 1
            N_Movement += 1;
            N_Displacement += 1;
            Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted = Movement(Patch_Radius, h, L, Beta, z_Displacement, Displacement, σ_p, λ_p, σ_w, λ_w, Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted, R_Cut, x, y, z)
        end
        if RN == 2
            N_Insertion += 1;
            #Pc_Grid, Pc_Grid_Sum, Pc_Grid_N = Grid_Excluded_Volume(2σ_p, h, L, Pc_Grid, Pc_Grid_Sum, Pc_Grid_N, x, y, z)
            Pc, Pc_Sum, Pc_N, x_Insertion, y_Insertion, z_Insertion = Random_Excluded_Volume(2σ_p, h, L, σ_w, Pc, Pc_Sum, Pc_N, x, y, z)
            #Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N = Analytic_Excluded_Volume(2σ_p, L, V, Pc_Analytic, Pc_Analytic_Sum, Pc_Analytic_N, x, y, z)
            if length(x_Insertion) > 0
                Energy, N_Insertion_Accepted, N_Insertion_Rejected = Insertion_Mezei(Patch_Radius, h, Beta, ChemPot, L, V, R_Cut, σ_p, λ_p, σ_w, λ_w, Energy, N_Insertion_Accepted, N_Insertion_Rejected, x, y, z, x_Insertion, y_Insertion, z_Insertion, Pc)
            else
                Energy, N_Insertion_Accepted, N_Insertion_Rejected = Insertion(Patch_Radius, h, L, V, ChemPot, Beta, R_Cut, σ_p, λ_p, σ_w, λ_w, x, y, z, Energy, N_Insertion_Accepted, N_Insertion_Rejected)
            end
        end
        if RN == 3 && length(x) > 1
            N_Removal += 1;
            if length(Pc) == 1
                Pc_Interpolation = Pc(collect(keys(Pc))[1])
            else
                if haskey(Pc, length(x) - 1)
                    Pc_Interpolation = Pc[length(x) - 1]
                else
                    Pc_Interpolation = Interpolation(Pc, length(x))
                end
            end
            if rand() > (1 - Pc_Interpolation)^1000
                Energy, N_Removal_Accepted, N_Removal_Rejected = Removal_Mezei(Patch_Radius, Pc_Interpolation, h, L, V, Beta, ChemPot, R_Cut, σ_p, λ_p, σ_w, λ_w, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z);
            else
                Energy, N_Removal_Accepted, N_Removal_Rejected = Removal(Patch_Radius, h, L, V, Beta, ChemPot, R_Cut, σ_p, λ_p, σ_w, λ_w, Energy, N_Removal_Accepted, N_Removal_Rejected, x, y, z)
            end
        end
        if i % MC_Measurement == 0
            if i > MC_Relaxation_Steps
                N_Measurements += 1;
                Energy_Array[N_Measurements] = Energy / length(x);
                println(Energy_File, "$N_Measurements\t$(round(Energy_Array[N_Measurements], digits = 6))")
                Energy_Sum += Energy / length(x);
                Average_Energy_Array[N_Measurements] = Energy_Sum / N_Measurements;
                println(Average_Energy_File, "$N_Measurements\t$(round(Average_Energy_Array[N_Measurements], digits = 6))")
                Density_Array[N_Measurements] = length(x) / V;
                println(Density_File, "$N_Measurements\t$(round(Density_Array[N_Measurements], digits = 6))")
                Density_Sum += length(x) / V;
                Average_Density_Array[N_Measurements] = Density_Sum / N_Measurements;
                println(Average_Density_File, "$N_Measurements\t$(round(Average_Density_Array[N_Measurements], digits = 6))")
                if N_Measurements > 1
                    σ_Energy[N_Measurements] = std(Energy_Array[1:N_Measurements]);
                    σ_Density[N_Measurements] = std(Density_Array[1:N_Measurements]);
                end
                g_x[N_Measurements, :] = Distribution(N_Bins, L, x);
                g_y[N_Measurements, :] = Distribution(N_Bins, L, y);
                g_z[N_Measurements, :] = Distribution(N_Bins, h + 2σ_w, z);
                PotentialFunction += Potential(Patch_Radius, N_Bins, L, h, σ_w, λ_w, x, y, z)
            end
            if N_Displacement_Accepted / N_Displacement > 0.55
                Displacement *= 1.05
                z_Displacement *= 1.05
            else
                Displacement *= 0.95
                z_Displacement *= 0.95
            end
            Displacement < 0.05 ? Displacement = 0.05 : nothing
            Displacement > L / 4. ? Displacement = L / 4. : nothing
            z_Displacement < 0.05 ? z_Displacement = 0.05 : nothing
            z_Displacement > (h + 2σ_w) / 8. ? z_Displacement = (h + 2σ_w) / 8. : nothing
            N_Displacement, N_Displacement_Accepted = 0, 0;
        end
    end
    close(Energy_File)
    close(Average_Energy_File)
    close(Density_File)
    close(Average_Density_File)
    close(Displacement_File)
    close(Insertion_File)
    close(Removal_File)
    
    Delta_xy = L / N_Bins;
    Delta_z = (h + 2σ_w) / N_Bins;
    r_xy = zeros(Float64, N_Bins)
    r_z = zeros(Float64, N_Bins)
    g_x *= (N_Bins / V);
    g_x_Normalized = g_x ./ Bulk_Density;
    g_x_Mean = mean(g_x, dims = 1)';
    g_x_Normalized_Mean = mean(g_x_Normalized, dims = 1)';
    g_x_Std = std(g_x, dims = 1)';
    g_x_Normalized_Std = std(g_x_Normalized, dims = 1)';
    g_x_Max = findmax(g_x_Mean)[1] + g_x_Std[findmax(g_x_Mean)[2][1]];
    g_x_Normalized_Max = findmax(g_x_Normalized_Mean)[1] + g_x_Normalized_Std[findmax(g_x_Normalized_Mean)[2][1]];
    g_y *= (N_Bins / V);
    g_y_Normalized = g_y ./ Bulk_Density;
    g_y_Mean = mean(g_y, dims = 1)';
    g_y_Normalized_Mean = mean(g_y_Normalized, dims = 1)';
    g_y_Std = std(g_y, dims = 1)';
    g_y_Normalized_Std = std(g_y_Normalized, dims = 1)';
    g_y_Max = findmax(g_y_Mean)[1] + g_y_Std[findmax(g_y_Mean)[2][1]];
    g_y_Normalized_Max = findmax(g_y_Normalized_Mean)[1] + g_y_Normalized_Std[findmax(g_y_Normalized_Mean)[2][1]];
    g_z *= (N_Bins / V);
    g_z_Normalized = g_z ./ Bulk_Density;
    g_z_Mean = mean(g_z, dims = 1)';
    g_z_Normalized_Mean = mean(g_z_Normalized, dims = 1)';
    g_z_Std = std(g_z, dims = 1)';
    g_z_Normalized_Std = std(g_z_Normalized, dims = 1)';
    g_z_Max = findmax(g_z_Mean)[1] + g_z_Std[findmax(g_z_Mean)[2][1]];
    g_z_Normalized_Max = findmax(g_z_Normalized_Mean)[1] + g_z_Normalized_Std[findmax(g_z_Normalized_Mean)[2][1]];
    PotentialFunction /= N_Measurements;

    g_x_File = open("$Output_Route/Distribution_x.dat", "w");
    println(g_x_File, "x\tDensity\tStandardDeviation\tNormalizedDensity\tStandardDeviation\n")
    g_y_File = open("$Output_Route/Distribution_y.dat", "w");
    println(g_y_File, "y\tDensity\tStandardDeviation\tNormalizedDensity\tStandardDeviation\n")
    g_z_File = open("$Output_Route/Distribution_z.dat", "w");
    println(g_z_File, "z\tDensity\tStandardDeviation\tNormalizedDensity\tStandardDeviation\n")
    Potential_File = open("$Output_Route/Potential_Function.dat", "w")
    println(Potential_File, "#z\t#U_wall(r)")

    @inbounds for i = 1:N_Bins
        r_xy[i] = round( - L / 2 + (i - 0.5) * Delta_xy, digits = 6);
        r_z[i] = round( - (h + 2σ_w) / 2 + (i - 0.5) * Delta_z, digits = 6);
        println(g_x_File, "$(r_xy[i])\t$(round(g_x_Mean[i], digits = 6))\t$(round(g_x_Std[i], digits = 6))\t$(round(g_x_Normalized_Mean[i], digits = 6))\t$(round(g_x_Normalized_Std[i], digits = 6))")
        println(g_y_File, "$(r_xy[i])\t$(round(g_y_Mean[i], digits = 6))\t$(round(g_y_Std[i], digits = 6))\t$(round(g_y_Normalized_Mean[i], digits = 6))\t$(round(g_y_Normalized_Std[i], digits = 6))")
        println(g_z_File, "$(r_z[i])\t$(round(g_z_Mean[i], digits = 6))\t$(round(g_z_Std[i], digits = 6))\t$(round(g_z_Normalized_Mean[i], digits = 6))\t$(round(g_z_Normalized_Std[i], digits = 6))")
        println(Potential_File, "$(r_z[i])\t$(round(PotentialFunction[i], digits = 6))")
    end
    close(g_x_File)
    close(g_y_File)
    close(g_z_File)
    close(Potential_File)

    g_x_Plot = plot(r_xy, g_x_Mean, ribbon = g_x_Std, ylim = (0, g_z_Max), legend = false, fillalpha = 0.2, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, title = "Slit Separation: $h", titlefontsize = 25, width = 3)
    g_x_Normalized_Plot = plot(r_xy, g_x_Normalized_Mean, ribbon = g_x_Normalized_Std, ylim = (0, g_z_Normalized_Max), legend = false, fillalpha = 0.2, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, title = "Slit Separation: $h", titlefontsize = 25, width = 3)
    hline!([1.], width = 3, linecolor = :black, linestyle = :dash)
    g_x_Plot_Individual = plot(r_xy, g_x_Mean, ribbon = g_x_Std, ylim = (0, g_z_Max), legend = false, fillalpha = 0.2, guidefontsize = 25, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "X Axis", ylabel = "Density", width = 3, size = [1920, 1080], widen = true, dpi = 300)
    #savefig(g_x_Plot_Individual, "$Output_Route/Density_x")
    g_x_Normalized_Plot_Individual = plot(r_xy, g_x_Normalized_Mean, ribbon = g_x_Normalized_Std, ylim = (0, g_z_Normalized_Max), legend = false, fillalpha = 0.2, guidefontsize = 25, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "X Axis", ylabel = "Density", width = 3, size = [1920, 1080], widen = true, dpi = 300)
    hline!([1.], width = 3, linecolor = :black, linestyle = :dash)
    #savefig(g_x_Normalized_Plot_Individual, "$Output_Route/NormalizedDensity_x")

    g_y_Plot = plot(r_xy, g_y_Mean, ribbon = g_y_Std, ylim = (0, g_z_Max), legend = false, fillalpha = 0.2, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, title = "Slit Separation: $h", titlefontsize = 25, width = 3)
    g_y_Normalized_Plot = plot(r_xy, g_y_Normalized_Mean, ribbon = g_y_Normalized_Std, ylim = (0, g_z_Normalized_Max), legend = false, fillalpha = 0.2, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, title = "Slit Separation: $h", titlefontsize = 25, width = 3)
    hline!([1.], width = 3, linecolor = :black, linestyle = :dash)
    g_y_Plot_Individual = plot(r_xy, g_y_Mean, ribbon = g_y_Std, ylim = (0, g_z_Max), legend = false, fillalpha = 0.2, guidefontsize = 25, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "Y Axis", title = "Slit Separation: $h", titlefontsize = 30, width = 3, size = [1920, 1080], widen = true, dpi = 300)
    #savefig(g_y_Plot_Individual, "$Output_Route/Density_y")
    g_y_Normalized_Plot_Individual = plot(r_xy, g_y_Normalized_Mean, ribbon = g_y_Normalized_Std, ylim = (0, g_z_Normalized_Max), legend = false, fillalpha = 0.2, guidefontsize = 25, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "Y Axis", title = "Slit Separation: $h", titlefontsize = 30, width = 3, size = [1920, 1080], widen = true, dpi = 300)
    hline!([1.], width = 3, linecolor = :black, linestyle = :dash)
    #savefig(g_y_Normalized_Plot_Individual, "$Output_Route/NormalizedDensity_y")

    g_z_Plot = plot(r_z, g_z_Mean, ribbon = g_z_Std, ylim = (0, g_z_Max), legend = false, fillalpha = 0.2, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, title = "Slit Separation: $h", titlefontsize = 25, width = 3)
    g_z_Normalized_Plot = plot(r_z, g_z_Normalized_Mean, ribbon = g_z_Normalized_Std, ylim = (0, g_z_Normalized_Max), legend = false, fillalpha = 0.2, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, title = "Slit Separation: $h", titlefontsize = 25, width = 3)
    hline!([1.], width = 3, linecolor = :black, linestyle = :dash)
    g_z_Plot_Individual = plot(r_z, g_z_Mean, ribbon = g_z_Std, ylim = (0, g_z_Max), legend = false, fillalpha = 0.2, guidefontsize = 25, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "Z Axis", width = 3, size = [1920, 1080], widen = true, dpi = 300)
    #savefig(g_z_Plot_Individual, "$Output_Route/Density_z")
    g_z_Normalized_Plot_Individual = plot(r_z, g_z_Normalized_Mean, ribbon = g_z_Normalized_Std, ylim = (0, g_z_Normalized_Max), legend = false, fillalpha = 0.2, guidefontsize = 25, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "Z Axis", width = 3, size = [1920, 1080], widen = true, dpi = 300)
    hline!([1.], width = 3, linecolor = :black, linestyle = :dash)
    #savefig(g_z_Normalized_Plot_Individual, "$Output_Route/NormalizedDensity_z")

    Density_Plots = plot(g_x_Plot_Individual, g_y_Plot_Individual, g_z_Plot_Individual, layout = (1, 3), link = :y, widen = true, size = [1920, 1080], dpi = 300)
    savefig(Density_Plots, "$Output_Route/Density_Profiles")
    Normalized_Density_Plots = plot(g_x_Normalized_Plot_Individual, g_y_Normalized_Plot_Individual, g_z_Normalized_Plot_Individual, layout = (1, 3), link = :y, widen = true, size = [1920, 1080], dpi = 300)
    savefig(Normalized_Density_Plots, "$Output_Route/Normalized_Density_Profiles")

    Potential_Plot = plot(r_z, PotentialFunction, legend = false, guidefontsize = 20, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "Z Axis", ylabel = "U_Wall(r)", title = "Slit Separation: $h", titlefontsize = 25, width = 3, size = [1920, 1080], widen = true, dpi = 300)
    savefig(Potential_Plot, "$Output_Route/Potential_Function")

    Pc_File = open("$Output_Route/Pc.dat", "w");
    println(Pc_File, "Density\tPc_Random\tPc_Grid\tPc_Analytic")
    Pc_Array = zeros(Float64, length(Pc) - 1)
    #Pc_Grid_Array = zeros(Float64, length(Pc) - 1)
    #Pc_Analytic_Array = zeros(Float64, length(Pc) - 1)
    @inbounds for i in keys(Pc)
        if i != 0
            Pc_Array[i] = Pc[i];
            #Pc_Grid_Array[i] = Pc_Grid[i];
            #Pc_Analytic_Array[i] = Pc_Analytic[i];
        end
    end
    Pc_Density = zeros(Float64, length(Pc_Array))
    for i = 1:length(Pc_Array)
        Pc_Density[i] = i / V;
        println(Pc_File, "$(Pc_Density[i])\t$(round(Pc_Array[i], digits = 6))")
        #println(Pc_File, "$i\t$(round(Pc_Array[i], digits = 6))\t$(round(Pc_Grid_Array[i], digits = 6))\t$(round(Pc_Analytic_Array[i], digits = 6))")
    end
    close(Pc_File)
    Pc_Plot = plot(Pc_Density, Pc_Array, legend = false, guidefontsize = 20, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "Density", ylabel = "Cavity Probability", ylims = (0, 1), title = "Slit Separation: $h", titlefontsize = 25, width = 3, size = [1920, 1080], widen = true, dpi = 300)
    #plot!(Pc_Grid_Array, label = "Grid", width = 3)
    #plot!(Pc_Analytic_Array, label = "Analytic", width = 3)
    savefig(Pc_Plot, "$Output_Route/Pc")

    D = CSV.read("$Output_Route/DisplacementAcceptance.dat", delim = "\t")
    Displacement_Plot = plot(D.N, D.Acceptance, legend = false, guidefontsize = 20, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "Measurement", ylabel = "Displacement Acceptance Percentage", ylims = (0, 100), title = "Slit Separation: $h", titlefontsize = 25, width = 3, size = [1920, 1080], widen = true, dpi = 300)
    savefig(Displacement_Plot, "$Output_Route/Displacement")

    D = CSV.read("$Output_Route/InsertionAcceptance.dat", delim = "\t")
    Insertion_Plot = plot(D.N, D.Acceptance, legend = false, guidefontsize = 20, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "Measurement", ylabel = "Insertion Acceptance Percentage", ylims = (0, 100), title = "Slit Separation: $h", titlefontsize = 25, width = 3, size = [1920, 1080], widen = true, dpi = 300)
    savefig(Insertion_Plot, "$Output_Route/Insertion")

    D = CSV.read("$Output_Route/RemovalAcceptance.dat", delim = "\t")
    Removal_Plot = plot(D.N, D.Acceptance, legend = false, guidefontsize = 20, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, xlabel = "Measurement", ylabel = "Removal Acceptance Percentage", ylims = (0, 100), title = "Slit Separation: $h", titlefontsize = 25, width = 3, size = [1920, 1080], widen = true, dpi = 300)
    savefig(Removal_Plot, "$Output_Route/Removal")
    ################################################################################
    println("< E / N > = $(round(mean(Energy_Array), digits = 6)) ± $(round(std(Energy_Array), digits = 6))")
    
    Energy_Plot_Individual = plot(Energy_Array, title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, xlabel = "Measurements", ylabel = "Energy [Unitless]", width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([mean(Energy_Array)], color = :black, width = 2, linestyle = :dash)
    Energy_Plot = plot(Energy_Array, title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm)
    hline!([mean(Energy_Array)], color = :black, width = 2, linestyle = :dash)

    Energy_Histogram_Individual = histogram(Energy_Array, legend = false, framestyle = :box, xlabel = "Energy [Unitless]", ylabel = "Normalized Frequency",  bins = 20, normalize = true, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    vline!([mean(Energy_Array)], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(mean(Energy_Array), std(Energy_Array)), width = 3, linecolor = :black)
    Energy_Histogram = histogram(Energy_Array, title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box,  bins = 20, normalize = true, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true)
    vline!([mean(Energy_Array)], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(mean(Energy_Array), std(Energy_Array)), width = 3, linecolor = :black)

    Average_Energy_Plot_Individual = plot(Average_Energy_Array, ribbon = σ_Energy, fillalpha = 0.2, legend = false, framestyle = :box, xlabel = "Measurements", ylabel = "Average Energy [Unitless]", width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([mean(Energy_Array)], color = :black, width = 2, linestyle = :dash)
    Average_Energy_Plot = plot(Average_Energy_Array, ribbon = σ_Energy, fillalpha = 0.2, title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([mean(Energy_Array)], color = :black, width = 2, linestyle = :dash)

    Energy_Plots = plot(Energy_Plot_Individual, Energy_Histogram_Individual, Average_Energy_Plot_Individual, layout = (@layout [a{0.3h} ; b c]))
    savefig(Energy_Plots, "$Output_Route/Energy")
    #################################################################################
    println("< N > = $(round(V*mean(Density_Array[1:end - 1]), digits = 6)) ± $(round(V*std(Density_Array[1:end - 1]), digits = 6))")
    println("< Density > = $(round(mean(Density_Array[1:end - 1]), digits = 6)) ± $(round(std(Density_Array[1:end - 1]), digits = 6))")
    
    Density_Plot_Individual = plot(Density_Array, title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, xlabel = "Measurements", ylabel = "Density [Unitless]", width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([mean(Density_Array)], color = :black, width = 2, linestyle = :dash)
    Density_Plot = plot(Density_Array, title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm)
    hline!([mean(Density_Array)], color = :black, width = 2, linestyle = :dash)
    
    Density_Histogram_Individual = histogram(Density_Array, legend = false, framestyle = :box, xlabel = "Density [Unitless]", ylabel = "Normalized Frequency",  bins = 20, normalize = true, width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    vline!([mean(Density_Array)], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(mean(Density_Array), std(Density_Array)), width = 3, linecolor = :black)
    Density_Histogram = histogram(Density_Array, title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box,  bins = 20, normalize = true, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true)
    vline!([mean(Density_Array)], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(mean(Density_Array), std(Density_Array)), width = 3, linecolor = :black)

    Average_Density_Plot_Individual = plot(Average_Density_Array, ribbon = σ_Density, fillalpha = 0.2, legend = false, framestyle = :box, xlabel = "Measurements", ylabel = "Average Density [Unitless]", width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([mean(Density_Array)], color = :black, width = 2, linestyle = :dash)
    Average_Density_Plot = plot(Average_Density_Array, ribbon = σ_Density, fillalpha = 0.2, title = "Slit Separation = $h", titlefontsize = 25, legend = false, framestyle = :box, width = 3, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([mean(Density_Array)], color = :black, width = 2, linestyle = :dash)

    Density_Plots = plot(Density_Plot_Individual, Density_Histogram_Individual, Average_Density_Plot_Individual, layout = (@layout [a{0.3h} ; b c]))
    savefig(Density_Plots, "$Output_Route/Density")
    #################################################################################
    Summary_File = open("$Output_Route/Summary.dat", "w")
    println(Summary_File, "$MC_Relaxation_Steps Relaxation Steps.")
    println(Summary_File, "$MC_Equilibrium_Steps Equilibrium Steps")
    println(Summary_File, "Measurements Every $MC_Measurement Steps.")
    println(Summary_File, "$N_Measurements Total Measurements.\n")
    println(Summary_File, "< E / N > = $(round(mean(Energy_Array), digits = 6)) ± $(round(std(Energy_Array), digits = 6))")
    println(Summary_File, "< N > = $(round(V*mean(Density_Array), digits = 6)) ± $(round(V*std(Density_Array), digits = 6))")
    println(Summary_File, "< Density > = $(round(mean(Density_Array), digits = 6)) ± $(round(std(Density_Array), digits = 6))")
    close(Summary_File)

    Povray_ini(h, L, ChemPot, T, λ_w, N_Image - 1, Patch_Percentage)
    Povray_Pov(h, L, ChemPot, T, σ_w, λ_w, Patch_Radius, Patch_Percentage)
    Povray_ini_Z_Axis(h, L, ChemPot, T, λ_w, N_Image - 1, Patch_Percentage)
    Povray_Pov_Z_Axis(h, L, ChemPot, T, σ_w, λ_w, Patch_Radius, Patch_Percentage)

    return mean(Energy_Array), std(Energy_Array), mean(Density_Array), std(Density_Array), Energy_Plot, Energy_Histogram, Average_Energy_Plot, Density_Plot, Density_Histogram, Average_Density_Plot, g_x_Plot, g_y_Plot, g_z_Plot, g_x_Normalized_Plot, g_y_Normalized_Plot, g_z_Normalized_Plot
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
    Delta_E = Energy_New - Energy_Old;
    if rand() < exp(-Beta * Delta_E)
        N_Movement_Accepted += 1;
        N_Displacement_Accepted += 1;
        Energy += Delta_E;
    else
        N_Movement_Rejected += 1;
        x[j] = x_Old;
        y[j] = y_Old;
        z[j] = z_Old;
    end
    return Energy, N_Movement_Accepted, N_Movement_Rejected, N_Displacement_Accepted
end

function Insertion(Patch_Radius::Float64, h::Float64, L::Float64, V::Float64, ChemPot::Float64, Beta::Float64, R_Cut::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, Energy::Float64, N_Insertion_Accepted::Int64, N_Insertion_Rejected::Int64)
    x_Insertion = L * (rand() - 0.5)
    y_Insertion = L * (rand() - 0.5)
    z_Insertion = (h + 2σ_w) * (rand() - 0.5)
    Energy_Insertion = Energy_Calculation(Patch_Radius, h, L, R_Cut, σ_p, λ_p, σ_w, λ_w, x_Insertion, y_Insertion, z_Insertion, x, y, z)
    if rand() < exp( Beta * (ChemPot - Energy_Insertion) + log(V / (length(x) + 1)) )
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion)
        append!(y, y_Insertion)
        append!(z, z_Insertion)
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
        append!(x, x_Insertion[j])
        append!(y, y_Insertion[j])
        append!(z, z_Insertion[j])
        Energy += Energy_Insertion;
    else
        N_Insertion_Rejected += 1;
    end
    return Energy, N_Insertion_Accepted, N_Insertion_Rejected
end

function Removal(Patch_Radius::Float64, h::Float64, L::Float64, V::Float64, Beta::Float64, ChemPot::Float64, R_Cut::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, Energy::Float64, N_Removal_Accepted::Int64, N_Removal_Rejected::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    j = rand(1:length(x));
    Energy_Removal = Energy_Calculation(Patch_Radius, h, L, R_Cut, σ_p, λ_p, σ_w, λ_w, x[j], y[j], z[j], x, y, z)
    if rand() < exp( Beta * (Energy_Removal - ChemPot) + log(length(x) / V) )
        N_Removal_Accepted += 1;
        deleteat!(x, j)
        deleteat!(y, j)
        deleteat!(z, j)
        Energy -= Energy_Removal;
    else
        N_Removal_Rejected += 1;
    end
    return Energy, N_Removal_Accepted, N_Removal_Rejected 
end

function Removal_Mezei(Patch_Radius::Float64, Pc_Interpolation::Float64, h::Float64, L::Float64, V::Float64, Beta::Float64, ChemPot::Float64, R_Cut::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, Energy::Float64, N_Removal_Accepted::Int64, N_Removal_Rejected::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    j = rand(1:length(x))
    Energy_Removal = Energy_Calculation(Patch_Radius, h, L, R_Cut, σ_p, λ_p, σ_w, λ_w, x[j], y[j], z[j], x, y, z)
    if rand() < ( length(x) / (V * Pc_Interpolation) ) * exp(Beta * (Energy_Removal - ChemPot))
        N_Removal_Accepted += 1;
        deleteat!(x, j)
        deleteat!(y, j)
        deleteat!(z, j)
        Energy -= Energy_Removal;
    else
        N_Removal_Rejected += 1;
    end
    return Energy, N_Removal_Accepted, N_Removal_Rejected
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

function Energy_Calculation(Patch_Radius::Float64, h::Float64, L::Float64, R_Cut::Float64, σ_p::Float64, λ_p::Float64, σ_w::Float64, λ_w::Float64, rx::Float64, ry::Float64, rz::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1})
    Energy = 0;
    for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
        Delta_x = rx - (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w);
        Delta_y = ry - (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w);
        Delta_z = abs(rz) - (h / 2 + σ_w)
        r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
        r_center = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)^2 + (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)^2)
        if r_center < Patch_Radius
            Energy += u_SquareWell(r2, σ_w, λ_w, 1.0);
        else
            Energy += u_SquareWell(r2, σ_w, λ_w, 0.0);
        end

        Delta_x = rx - (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w);
        Delta_y = ry - (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w);
        Delta_z = abs(rz) - (h / 2 + σ_w)
        r2 = Delta_x^2 + Delta_y^2 + Delta_z^2;
        r_center = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)^2 + ((- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w))^2 );
        if r_center < Patch_Radius
            Energy += u_SquareWell(r2, σ_w, λ_w, 1.0);
        else
            Energy += u_SquareWell(r2, σ_w, λ_w, 0.0);
        end
        if Energy == Inf
            return Energy
        end
    end
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
    end
    return Energy
end

function PeriodicBoundaryConditions(L::Float64, x::Float64)
    if x < - L / 2.
        x += L;
    elseif x > L / 2.
        x -= L;
    end
    return x
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

    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/h_$(h)/Positions"
    Pov_File = open("$Output_Route/Pore_Animation.pov", "w");
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
end

function Povray_Pov_Z_Axis(h::Float64, L::Float64, ChemPot::Float64, T::Float64, σ_w::Float64, λ_w::Float64, Patch_Radius::Float64, Patch_Percentage::Int64)
    Particle_r, Particle_g, Particle_b, Particle_t = 174/255, 214/255, 241/255, 0;
    Slab_Attractive_r, Slab_Attractive_g, Slab_Attractive_b, Slab_Attractive_t = 241/255, 196/255, 15/255, 0; 
    Slab_Neutral_r, Slab_Neutral_g, Slab_Neutral_b, Slab_Neutral_t  = 142/255, 68/255, 173/255, 0; 
    Slab_Repulsive_r, Slab_Repulsive_g, Slab_Repulsive_b, Slab_Repulsive_t = 183/255, 28/255, 28/255, 0;
    Wall_r, Wall_g, Wall_b, Wall_t = 174/255, 214/255, 241/255, 0.5;
    
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/h_$(h)/Positions"
    Pov_File = open("$Output_Route/Pore_Z_Axis_Animation.pov", "w");
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
end

function Povray_ini(h::Float64, L::Float64, ChemPot::Float64, T::Float64, λ_w::Float64, Frames::Int64, Patch_Percentage::Int64)
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/h_$(h)/Positions"
    mkpath("$Output_Route")
    Ini_File = open("$Output_Route/Pore_Animation.ini", "w");
    println(Ini_File, "Input_File_Name = $Output_Route/Pore_Animation.pov")
    println(Ini_File, "Output_File_Name = $Output_Route/")
    println(Ini_File, "+W800 +H800\n")
    println(Ini_File, "Initial_Frame = 1")
    println(Ini_File, "Final_Frame = $Frames")
    println(Ini_File, "Initial_Clock = 1")
    println(Ini_File, "Final_Clock = $Frames\n")
    println(Ini_File, "Cyclic_Animation = off")
    close(Ini_File)
end

function Povray_ini_Z_Axis(h::Float64, L::Float64, ChemPot::Float64, T::Float64, λ_w::Float64, Frames::Int64, Patch_Percentage::Int64)
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/h_$(h)/Positions"
    mkpath("$Output_Route")
    Ini_File = open("$Output_Route/Pore_Z_Axis_Animation.ini", "w");
    println(Ini_File, "Input_File_Name = $Output_Route/Pore_Z_Axis_Animation.pov")
    println(Ini_File, "Output_File_Name = $Output_Route/")
    println(Ini_File, "+W800 +H800\n")
    println(Ini_File, "Initial_Frame = 1")
    println(Ini_File, "Final_Frame = $Frames")
    println(Ini_File, "Initial_Clock = 1")
    println(Ini_File, "Final_Clock = $Frames\n")
    println(Ini_File, "Cyclic_Animation = off")
    close(Ini_File)
end

function Povray_Slab(L::Float64, ChemPot::Float64, T::Float64, σ_w::Float64, λ_w::Float64, Patch_Radius::Float64, Patch_Percentage::Int64)
    Slab_Attractive_r, Slab_Attractive_g, Slab_Attractive_b, Slab_Attractive_t = 241/255, 196/255, 15/255, 0; 
    Slab_Neutral_r, Slab_Neutral_g, Slab_Neutral_b, Slab_Neutral_t  = 142/255, 68/255, 173/255, 0; 
    Slab_Repulsive_r, Slab_Repulsive_g, Slab_Repulsive_b, Slab_Repulsive_t = 183/255, 28/255, 28/255, 0;
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/"
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
    #if λ_w > 1.
    #    println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n\t\n#end")
    #else
    #    println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Neutral_r, $Slab_Neutral_g, $Slab_Neutral_b, $Slab_Neutral_t> }\n\tno_shadow}\n\t\n#end")
    #end
    println(Pov_File, """#fopen File_Slab "$Output_Route/Slab.xyz" read""")
    println(Pov_File, "#while (defined (File_Slab)) \n\t#read (File_Slab, rx, ry, rz)\n\tSlab(rx, ry, 0, $σ_w, $λ_w)\n#end\n#fclose File_Slab")
    close(Pov_File)
end

function Cycled_Mezei()
    ChemPot = -4.150918;
    L = 20.;
    T = 1.5;
    Bulk_Density = 0.1;
    H = range(2., 10., step = 1.);
    σ_p, λ_p = 0.5, 1.5;
    σ_w, λ_w = 0.5, 1.5;

    Patch_Percentage = 0;

    Mean_Energy = zeros(Float64, length(H));
    Std_Energy = zeros(Float64, length(H));
    Mean_Density = zeros(Float64, length(H));
    Std_Density = zeros(Float64, length(H));
    j = 1;

    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/Patch_$Patch_Percentage%/"
    mkpath("$Output_Route")
    p_x = Array{Any, 1}(undef, length(H))
    p_y = Array{Any, 1}(undef, length(H))
    p_z = Array{Any, 1}(undef, length(H))
    p_Normalized_x = Array{Any, 1}(undef, length(H))
    p_Normalized_y = Array{Any, 1}(undef, length(H))
    p_Normalized_z = Array{Any, 1}(undef, length(H))
    p_Energy = Array{Any, 1}(undef, length(H))
    p_Histogram_Energy = Array{Any, 1}(undef, length(H))
    p_Average_Energy = Array{Any, 1}(undef, length(H))
    p_Density = Array{Any, 1}(undef, length(H))
    p_Histogram_Density = Array{Any, 1}(undef, length(H))
    p_Average_Density = Array{Any, 1}(undef, length(H))

    Density_File = open("$Output_Route/Density_T_$(T).dat", "w");
    println(Density_File, "h\tDensity\tStandardDeviation")
    Energy_File = open("$Output_Route/Energy_T_$(T).dat", "w");
    println(Energy_File, "h\tEnergy\tStandardDeviation")
    for h in H
        Mean_Energy[j], Std_Energy[j], Mean_Density[j], Std_Density[j], p_Energy[j], p_Histogram_Energy[j], p_Average_Energy[j], p_Density[j], p_Histogram_Density[j], p_Average_Density[j], p_x[j], p_y[j], p_z[j], p_Normalized_x[j], p_Normalized_y[j], p_Normalized_z[j] = Mezei(ChemPot, h, L, T, σ_p, λ_p, σ_w, λ_w, j, length(H), Patch_Percentage, Bulk_Density);
        println(Density_File, "$h\t$(Mean_Density[j])\t$(Std_Density[j])")
        println(Energy_File, "$h\t$(Mean_Energy[j])\t$(Std_Energy[j])")
        j += 1;
    end
    close(Density_File)
    close(Energy_File)

    y = ones(3);
    title = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, y[2], Plots.text("T = $T, \\mu = $(round(ChemPot, digits = 3)), Bulk Density = $Bulk_Density", :black, 50)), axis=false, grid = false, leg=false)
    y_axis = plot(guide_position = :right, left_margin = -30mm, right_margin = 20mm, ylabel = "Density", guidefontsize = 30, grid = false, axis = false)

    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("X Axis", :black, 30)), axis=false, grid = false, leg=false)
    x_Plot = plot(y_axis, title, p_x[1], p_x[2], p_x[3], p_x[4], p_x[5], p_x[6], p_x[7], p_x[8], p_x[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(x_Plot, "$Output_Route/Distribution_x_$T")
    
    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("Y Axis", :black, 30)), axis=false, grid = false, leg=false)
    y_Plot = plot(y_axis, title, p_y[1], p_y[2], p_y[3], p_y[4], p_y[5], p_y[6], p_y[7], p_y[8], p_y[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(y_Plot, "$Output_Route/Distribution_y_$T")

    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("Z Axis", :black, 30)), axis=false, grid = false, leg=false)
    z_Plot = plot(y_axis, title, p_z[1], p_z[2], p_z[3], p_z[4], p_z[5], p_z[6], p_z[7], p_z[8], p_z[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(z_Plot, "$Output_Route/Distribution_z_$T")

    y_axis = plot(guide_position = :right, left_margin = -30mm, right_margin = 20mm, ylabel = "Normalized Density", guidefontsize = 30, grid = false, axis = false)
    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("X Axis", :black, 30)), axis=false, grid = false, leg=false)
    Normalized_x_Plot = plot(y_axis, title, p_Normalized_x[1], p_Normalized_x[2], p_Normalized_x[3], p_Normalized_x[4], p_Normalized_x[5], p_Normalized_x[6], p_Normalized_x[7], p_Normalized_x[8], p_Normalized_x[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Normalized_x_Plot, "$Output_Route/NormalizedDistribution_x_$T")
    
    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("Y Axis", :black, 30)), axis=false, grid = false, leg=false)
    Normalized_y_Plot = plot(y_axis, title, p_Normalized_y[1], p_Normalized_y[2], p_Normalized_y[3], p_Normalized_y[4], p_Normalized_y[5], p_Normalized_y[6], p_Normalized_y[7], p_Normalized_y[8], p_Normalized_y[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Normalized_y_Plot, "$Output_Route/NormalizedDistribution_y_$T")

    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("Z Axis", :black, 30)), axis=false, grid = false, leg=false)
    Normalized_z_Plot = plot(y_axis, title, p_Normalized_z[1], p_Normalized_z[2], p_Normalized_z[3], p_Normalized_z[4], p_Normalized_z[5], p_Normalized_z[6], p_Normalized_z[7], p_Normalized_z[8], p_Normalized_z[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Normalized_z_Plot, "$Output_Route/NormalizedDistribution_z_$T")

    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("Measurements", :black, 30)), axis=false, grid = false, leg=false)
    y_axis = plot(guide_position = :right, left_margin = -30mm, right_margin = 20mm, ylabel = "Energy [Unitless]", guidefontsize = 30, grid = false, axis = false)
    Energy_Plot = plot(y_axis, title, p_Energy[1], p_Energy[2], p_Energy[3], p_Energy[4], p_Energy[5], p_Energy[6], p_Energy[7], p_Energy[8], p_Energy[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Energy_Plot, "$Output_Route/Energy_$T")

    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("Energy [Unitless]", :black, 30)), axis=false, grid = false, leg=false)
    y_axis = plot(guide_position = :right, left_margin = -30mm, right_margin = 20mm, ylabel = "Normalized Frequency", guidefontsize = 30, grid = false, axis = false)
    Energy_Histogram_Plot = plot(y_axis, title, p_Histogram_Energy[1], p_Histogram_Energy[2], p_Histogram_Energy[3], p_Histogram_Energy[4], p_Histogram_Energy[5], p_Histogram_Energy[6], p_Histogram_Energy[7], p_Histogram_Energy[8], p_Histogram_Energy[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Energy_Histogram_Plot, "$Output_Route/Energy_Histogram_$T")

    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("Measurements", :black, 30)), axis=false, grid = false, leg=false)
    y_axis = plot(guide_position = :right, left_margin = -30mm, right_margin = 20mm, ylabel = "Average Energy [Unitless]", guidefontsize = 30, grid = false, axis = false)
    Average_Energy_Plot = plot(y_axis, title, p_Average_Energy[1], p_Average_Energy[2], p_Average_Energy[3], p_Average_Energy[4], p_Average_Energy[5], p_Average_Energy[6], p_Average_Energy[7], p_Average_Energy[8], p_Average_Energy[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Average_Energy_Plot, "$Output_Route/Average_Energy_$T")

    y_axis = plot(guide_position = :right, left_margin = -30mm, right_margin = 20mm, ylabel = "Density", guidefontsize = 30, grid = false, axis = false)
    Density_Plot = plot(y_axis, title, p_Density[1], p_Density[2], p_Density[3], p_Density[4], p_Density[5], p_Density[6], p_Density[7], p_Density[8], p_Density[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Density_Plot, "$Output_Route/Density_$T")

    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("Density", :black, 30)), axis=false, grid = false, leg=false)
    y_axis = plot(guide_position = :right, left_margin = -30mm, right_margin = 20mm, ylabel = "Normalized Frequency", guidefontsize = 30, grid = false, axis = false)
    Density_Histogram_Plot = plot(y_axis, title, p_Histogram_Density[1], p_Histogram_Density[2], p_Histogram_Density[3], p_Histogram_Density[4], p_Histogram_Density[5], p_Histogram_Density[6], p_Histogram_Density[7], p_Histogram_Density[8], p_Histogram_Density[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Density_Histogram_Plot, "$Output_Route/Density_Histogram_$T")

    x_axis = Plots.scatter(y, marker = 0, markeralpha = 0, annotations=(2, 1.5*y[2], Plots.text("Measurements", :black, 30)), axis=false, grid = false, leg=false)
    y_axis = plot(guide_position = :right, left_margin = -30mm, right_margin = 20mm, ylabel = "Average Density", guidefontsize = 30, grid = false, axis = false)
    Average_Density_Plot = plot(y_axis, title, p_Average_Density[1], p_Average_Density[2], p_Average_Density[3], p_Average_Density[4], p_Average_Density[5], p_Average_Density[6], p_Average_Density[7], p_Average_Density[8], p_Average_Density[9], x_axis, layout = (@layout [ a{0.005w} [b{0.05h} ; grid(3,3) ; c{0.05h} ]]), size = [1920, 1080], dpi = 300)
    savefig(Average_Density_Plot, "$Output_Route/Average_Density_$T")

    #Total_Average_Density_Plot = plot(H, Mean_Density, yerror = Std_Density, ylims = (0, findmax(Mean_Density)[1] + Std_Density[findmax(Mean_Density)[2]] ), markerstrokewidth = 3, title = "T = $T, \\mu = $(round(ChemPot, digits = 3)), Bulk Density = $Bulk_Density", titlefontsize = 25, legend = false, framestyle = :box, xlabel = "Pore Separation", ylabel = "Average Density [Unitless]", width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    #savefig(Total_Average_Density_Plot, "$Output_Route/AverageDensity")

    #Total_Average_Energy_Plot = plot(H, Mean_Energy, yerror = Std_Energy, markerstrokewidth = 3, title = "T = $T, \\mu = $(round(ChemPot, digits = 3)), \\rho_{Bulk} = $Bulk_Density", titlefontsize = 25, legend = false, framestyle = :box, xlabel = "Pore Separation", ylabel = "Average Energy [Unitless]", width = 3, guidefontsize = 20, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    #savefig(Total_Average_Energy_Plot, "$Output_Route/AverageEnergy")

    N_Slab = 2 * convert(Int64, ceil(L / (2σ_w) + 3 + 1) * ceil(L / (√3 * 2σ_w) + √3))
    File_Slabs_Avogadro = open("$Output_Route/Slab_Avogadro.xyz", "w+")
    println(File_Slabs_Avogadro, "$N_Slab\n")
    File_Slabs_Povray = open("$Output_Route/Slab.xyz", "w+")
    for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
        r = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)^2 + (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)^2)
        if r < 5.
            println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
        else
            println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
        end
        r = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)^2 + (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)^2)
        if r < 5.
            println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
        else
            println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
        end
            if i_x == ceil(L / (2σ_w) + 3 + 1) && i_y == ceil(L / (√3 * 2σ_w) + √3)
            println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w),\t0.,")
            println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w),\t0.")
        else
            println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w),\t0.,")
            println(File_Slabs_Povray, "$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w),\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w),\t0.,")
        end
    end
    close(File_Slabs_Avogadro)
    close(File_Slabs_Povray)

    Patch_Radius = round(sqrt(L^2 * Patch_Percentage / 314) , digits = 6)
    Povray_Slab(L, ChemPot, T, σ_w, λ_w, Patch_Radius, Patch_Percentage)

    for h in H
        run(`povray $Output_Route/h_$h/Positions/Pore_Animation.ini`)
        run(`povray $Output_Route/h_$h/Positions/Pore_Z_Axis_Animation.ini`)
    end
    run(`povray $Output_Route/Pore_Slab.pov +W800 +H800`)

end

@time Cycled_Mezei()