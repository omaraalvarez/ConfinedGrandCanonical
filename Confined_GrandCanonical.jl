using Statistics;
using Plots; gr();
using Plots.PlotMeasures;
using Test;
using LaTeXStrings;
using StatsPlots;
using Distributions;
using ArgParse;

function Parse_Commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "μ"
            help = "Chemical Potental";
            arg_type = Float64;
            required = true;
        "T"
            help = "Temperature";
            arg_type = Float64;
            required = true;
        "h"
            help = "Slit's Separation";
            arg_type = Float64;
            required = true;
        "--ρ_Bulk", "-D"
            arg_type = Float64;
            default = 1.;
            help = "Bulk fluid's density";
            required = false;
        "--Configurations", "-C"
            help = "Enable/Disable Saving Configuraitons";
            arg_type = Bool;
            default = false;
        "--PovRay", "-P"
            help = "Enable/Disable PovRay Animations";
            arg_type = Bool;
            default = false;
    end
    return parse_args(s)
end

function Confined_GrandCanonical_MonteCarlo(μ::Type, T::Type, h::Type, Bulk_Density::Type = 1., Configurations::Bool = false, PovRay::Bool = false, R_Cut::Type = 3.) where {Type <: Real}
    ####################################################     CONFIGURATIONAL STEPS   ################################################################
    println("\n\tCONFINED GRAND CANONICAL MONTE CARLO")
    MC_Relaxation_Steps = 20_000;
    MC_Equilibrium_Steps = 250_000;
    MC_Measurement = 10;
    MC_Steps = MC_Equilibrium_Steps + MC_Relaxation_Steps;
    ##################################### VARIABLE INITIALIZATION ###########################
    L, σ_w, λ_w, σ_p, λ_p, Delta_Bins = 20., 1.0, 1.5, 1.0, 1.5, 0.02;
    V, Beta, Nxy_Bins, Nz_Bins, Equilibrium = h * L^2, 1. / T, convert(Int64, ceil(L / Delta_Bins)), convert(Int64, ceil(h / Delta_Bins)), false;
    N_Id, N_Image = convert(Int64, round(0.5 * V)), 1;
    Pc, Pc_Sum, Pc_N = Dict{Int64, Float64}(), Dict{Int64, Float64}(), Dict{Int64, Int64}();
    Z_Displacement, Displacement, N_Displacement, N_Displacement_Accepted = h / 8., L / 8., 0, 0;
    N_Insertion, N_Insertion_Accepted = 0, 0;
    N_Removal, N_Removal_Accepted = 0, 0;
    Energy, Density, N_Measurements = 0., 0., 0;
    Energy_Array, Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    Mean_Energy_Array, Mean_Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    STD_Energy_Array, STD_Density_Array = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) ), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ) );
    g_x, g_y, g_z = zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ), Nxy_Bins), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ), Nxy_Bins), zeros(Float64, convert(Int64, ceil(MC_Equilibrium_Steps / MC_Measurement) ), Nz_Bins);
    ####################################### OUTPUT ROUTE, INITIAL POSITIONS AND INITIAL ENERGY ###########################
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(μ, digits = 3))/h_$(round(h, digits = 2))"
    mkpath("$Output_Route")
    PovRay ? mkpath("$Output_Route/Positions") : nothing
    Configurations ? mkpath("$Output_Route/Configurations") : nothing
    x, y, z = Float64[], Float64[], Float64[];
    ####################################################     SLAB CREATION   ################################################################
    if PovRay
        File_Slabs_Povray = open("$Output_Route/Positions/Slab.xyz", "w+")
        for i_y = 1:ceil(L / (√3 * σ_w) + √3), i_x = 1:ceil(L / σ_w + 3 + 1)
            if i_x == ceil(L / σ_w + 3 + 1) && i_y == ceil(L / (√3 * σ_w) + √3)
                println(File_Slabs_Povray, "$(- L / 2 - 1.5σ_w + (i_x - 1) * σ_w),\t$(- L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w),\t$(h / 2),")
                println(File_Slabs_Povray, "$(- L / 2 - 1.5σ_w + (i_x - 1) * σ_w + σ_w/2),\t$(- L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w + √3/2 * σ_w),\t$(h / 2)")
            else
                println(File_Slabs_Povray, "$(- L / 2 - 1.5σ_w + (i_x - 1) * σ_w),\t$(- L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w),\t$(h / 2),")
                println(File_Slabs_Povray, "$(- L / 2 - 1.5σ_w + (i_x - 1) * σ_w + σ_w/2),\t$(- L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w + √3/2 * σ_w),\t$(h / 2),")
            end
        end
        close(File_Slabs_Povray)
    end
    ####################################################    SIMULATION LOOP   ################################################################
    @inbounds for k = 1:MC_Steps
        @test all(Array(abs.(x)) .<= L / 2.)
        @test all(Array(abs.(y)) .<= L / 2.)
        @test all(Array(abs.(z)) .<= h / 2.)
        ####################################### PRINTS SIMULATION PROGRESS TO SCREEN ########################################
        if k < MC_Relaxation_Steps && k % .01MC_Relaxation_Steps == 0
            println("$(convert(Int64, 100k / MC_Relaxation_Steps))% Relaxation: [μ = $μ, L = $L, h = $h, T = $T]")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("N = $(length(x)) Particles")
            println("ρ = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))\tMax Displacement Z Axis = $(round(Z_Displacement, digits = 6))")
            println("Movements: $N_Displacement")
            println("   Accepted: $N_Displacement_Accepted ($(round(100N_Displacement_Accepted/N_Displacement, digits = 2))%)\tRejected: $(N_Displacement - N_Displacement_Accepted) ($(round(100(N_Displacement - N_Displacement_Accepted)/N_Displacement, digits = 2))%)")
            println("Insertions: $N_Insertion")
            println("   Accepted: $N_Insertion_Accepted ($(round(100N_Insertion_Accepted / N_Insertion, digits = 2))%)\tRejected: $(N_Insertion - N_Insertion_Accepted) ($(round(100(N_Insertion - N_Insertion_Accepted) / N_Insertion, digits = 2))%)")
            println("Removal: $N_Removal")
            println("   Accepted: $N_Removal_Accepted ($(round(100N_Removal_Accepted / N_Removal, digits = 2))%)\tRejected: $(N_Removal - N_Removal_Accepted) ($(round(100(N_Removal - N_Removal_Accepted) / N_Removal, digits = 2))%)\n")
            N_Displacement, N_Displacement_Accepted = 0, 0;
            N_Insertion, N_Insertion_Accepted = 0, 0;
            N_Removal, N_Removal_Accepted = 0, 0;
        end

        if k > MC_Relaxation_Steps && k % .01MC_Equilibrium_Steps == 0
            println("$(convert(Int64, ceil(100(k - MC_Relaxation_Steps) / MC_Equilibrium_Steps)))% Equilibrium: [μ = $μ, L = $L, h = $h, T = $T, $(N_Measurements + 1) Measurements]")
            println("U / N = $(round(Energy / length(x), digits = 6))")
            println("N = $(length(x)) Particles")
            println("ρ = $(round(length(x) / V, digits = 6))")
            println("Max Displacement = $(round(Displacement, digits = 6))\tMax Displacement Z Axis = $(round(Z_Displacement, digits = 6))")
            println("Movements: $N_Displacement")
            println("   Accepted: $N_Displacement_Accepted ($(round(100N_Displacement_Accepted/N_Displacement, digits = 2))%)\tRejected: $(N_Displacement - N_Displacement_Accepted) ($(round(100(N_Displacement - N_Displacement_Accepted)/N_Displacement, digits = 2))%)")
            println("Insertions: $N_Insertion")
            println("   Accepted: $N_Insertion_Accepted ($(round(100N_Insertion_Accepted / N_Insertion, digits = 2))%)\tRejected: $(N_Insertion - N_Insertion_Accepted) ($(round(100(N_Insertion - N_Insertion_Accepted) / N_Insertion, digits = 2))%)")
            println("Removal: $N_Removal")
            println("   Accepted: $N_Removal_Accepted ($(round(100N_Removal_Accepted / N_Removal, digits = 2))%)\tRejected: $(N_Removal - N_Removal_Accepted) ($(round(100(N_Removal - N_Removal_Accepted) / N_Removal, digits = 2))%)\n")
            N_Displacement, N_Displacement_Accepted = 0, 0;
            N_Insertion, N_Insertion_Accepted = 0, 0;
            N_Removal, N_Removal_Accepted = 0, 0;
        end

        k == MC_Relaxation_Steps ? println("- FINISHED RELAXATION STEPS\n") : nothing
        k == MC_Relaxation_Steps ? Equilibrium = true : nothing
        ####################################################    SAVES POSITIONS OF THE MOLECULES FOR THE POVRAY IMAGES  ################################################################
        if PovRay
            if k > MC_Relaxation_Steps && k % .1MC_Equilibrium_Steps == 0
                Positions_File = open("$Output_Route/Positions/Pos_$N_Image.xyz", "w");
                @inbounds for i = 1:length(x)
                    i != length(x) ? println(Positions_File, "$(x[i]),\t$(y[i]),\t$(z[i]),") : println(Positions_File, "$(x[i]),\t$(y[i]),\t$(z[i])")
                end
                close(Positions_File)
                N_Image += 1
            end
        end

        @inbounds for i = 1:N_Id
        RN = rand(1:3);
            ####################################################    MOVEMENT   ################################################################
            if RN == 1 && length(x) > 1
                N_Displacement += 1;
                Energy, N_Displacement_Accepted = Movement(h, L, Beta, Z_Displacement, Displacement, Energy, N_Displacement_Accepted, x, y, z, σ_p, λ_p, σ_w, λ_w)
            end
            ####################################################    INSERTION   ################################################################
            if RN == 2
                N_Insertion += 1;
                Pc, Pc_Sum, Pc_N, x_Insertion, y_Insertion, z_Insertion = Random_Excluded_Volume(Equilibrium, h, L, Pc, Pc_Sum, Pc_N, x, y, z, σ_p, σ_w)
                if length(x_Insertion) > 0
                    Energy, N_Insertion_Accepted = Insertion_Mezei(h, L, μ, Beta, Energy, N_Insertion_Accepted, x, y, z, x_Insertion, y_Insertion, z_Insertion, Pc, σ_p, λ_p, σ_w, λ_w)
                else
                    Energy, N_Insertion_Accepted = Insertion(h, L, μ, Beta, Energy, N_Insertion_Accepted, x, y, z, σ_p, λ_p, σ_w, λ_w)
                end
            end
            ####################################################     REMOVAL    ################################################################
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
                if rand() > (1 - Pc_Interpolation)^(Equilibrium ? 200 : 1000)
                    Energy, N_Removal_Accepted = Removal_Mezei(h, L, μ, Beta, Energy, N_Removal_Accepted, x, y, z, Pc_Interpolation, σ_p, λ_p, σ_w, λ_w);
                else
                    Energy, N_Removal_Accepted = Removal(h, L, μ, Beta, Energy, N_Removal_Accepted, x, y, z, σ_p, λ_p, σ_w, λ_w)
                end
            end
        end

        ####################################################    MEASUREMENT SECTION   ################################################################
        if k % MC_Measurement == 0
            if k > MC_Relaxation_Steps
                N_Measurements += 1;
                Energy_Array[N_Measurements] = Energy / length(x);
                Density_Array[N_Measurements] = length(x) / V;
                Mean_Energy_Array[N_Measurements] = mean(Energy_Array[1:N_Measurements]);
                Mean_Density_Array[N_Measurements] = mean(Density_Array[1:N_Measurements])
                if N_Measurements > 1
                    STD_Energy_Array[N_Measurements] = std(Energy_Array[1:N_Measurements]);
                    STD_Density_Array[N_Measurements] = std(Density_Array[1:N_Measurements]);
                end
                g_x[N_Measurements, :] = Distribution(Nxy_Bins, L, x);
                g_y[N_Measurements, :] = Distribution(Nxy_Bins, L, y);
                g_z[N_Measurements, :] = Distribution(Nz_Bins, h, z);

                ####################################################    SAVES POSITIONS OF THE MOLECULES    ################################################################
                if Configurations
                    Configurations_File = open("$Output_Route/Configurations/Conf_$N_Measurements.xyz", "w");
                    println(Configurations_File, "x\ty\tz")
                    @inbounds for i = 1:length(x)
                        println(Configurations_File, "$(x[i])\t$(y[i])\t$(z[i])")
                    end
                    close(Configurations_File)
                end
            end

            ####################################################    MAX DISPLACEMENT CONTROL  ################################################################
            N_Displacement_Accepted / N_Displacement > 0.55 ? Displacement *= 1.05 : Displacement *= 0.95
            N_Displacement_Accepted / N_Displacement > 0.55 ? Z_Displacement *= 1.05 : Z_Displacement *= 0.95
            Displacement < 0.05 ? Displacement = 0.05 : (Displacement > L / 4. ? _Displacement = L / 4. : nothing)
            Z_Displacement < 0.05 ? Z_Displacement = 0.05 : (Z_Displacement > h / 8. ? Z_Displacement = h / 8. : nothing)
        end
    end
    ####################################################### END OF SIMULATION CYCLES #############################################

    ############################################################### SUMMARY FILE #################################################
    println("< U* / N > = $(round(Mean_Energy_Array[N_Measurements], digits = 6)) ± $(round(STD_Energy_Array[N_Measurements], digits = 6))")
    println("< ρ* > = $(round(Mean_Density_Array[N_Measurements], digits = 6)) ± $(round(STD_Density_Array[N_Measurements], digits = 6))")
    Summary_File = open("$Output_Route/Summary.dat", "w+")
    println(Summary_File, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   INPUT   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    println(Summary_File, " μ = $μ\nL = $L\th = $h\tV = $V\nT = $T\n$MC_Relaxation_Steps Relaxation Steps.\n$MC_Equilibrium_Steps Equilibrium Steps.\tMeasurements every $MC_Measurement steps.")
    println(Summary_File, "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   OUTPUT   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    println(Summary_File, "< U* / N > = $(round(Mean_Energy_Array[N_Measurements], digits = 6)) ± $(round(STD_Energy_Array[N_Measurements], digits = 6))")
    println(Summary_File, "< ρ* > = $(round(Mean_Density_Array[N_Measurements], digits = 6)) ± $(round(STD_Density_Array[N_Measurements], digits = 6))")
    close(Summary_File)
    ################################################################ OUTPUT ##############################################################
    Energy_File = open("$Output_Route/Energy.dat", "w+");
    Density_File = open("$Output_Route/Density.dat", "w+");
    println(Energy_File, "Step\tEnergy\tMean_Energy\tStd_Energy")
    println(Density_File, "Step\tDensity\tMean_Density\tStd_Density")
    for i = 1:N_Measurements
        println(Energy_File, "$i\t$(round(Energy_Array[i], digits = 6))\t$(round(Mean_Energy_Array[i], digits = 6))\t$(round(STD_Energy_Array[i], digits = 6))")
        println(Density_File, "$i\t$(round(Density_Array[i], digits = 6))\t$(round(Mean_Density_Array[i], digits = 6))\t$(round(STD_Density_Array[i], digits = 6))")
    end
    close(Energy_File), close(Density_File)    
    #################################################################   DISTRIBUTION PROFILES  ################################################################
    Delta_xy, Delta_z = L / Nxy_Bins, h / Nz_Bins;
    r_xy, r_z = zeros(Float64, Nxy_Bins), zeros(Float64, Nz_Bins);
    g_x *= (Nxy_Bins / (V * Bulk_Density));
    g_x_Mean, g_x_Std = mean(g_x, dims = 1)', std(g_x, dims = 1)';
    g_y *= (Nxy_Bins / (V * Bulk_Density));
    g_y_Mean, g_y_Std = mean(g_y, dims = 1)', std(g_y, dims = 1)';
    g_z *= (Nz_Bins / (V * Bulk_Density));
    g_z_Mean, g_z_Std = mean(g_z, dims = 1)', std(g_z, dims = 1)';
    g_z_Max = findmax(g_z_Mean)[1] + g_z_Std[findmax(g_z_Mean)[2][1]];

    Density_Distribution_File = open("$Output_Route/Density_Distribution.dat", "w");
    Density_Distribution_Z_File = open("$Output_Route/Density_Distribution_Z.dat", "w");
    println(Density_Distribution_File, "r_xy\tDensity_x\tSTD_x\tDensity_y\tSTD_y")
    println(Density_Distribution_Z_File, "r_z\tDensity_z\tSTD_z")
    @inbounds for i = 1:Nxy_Bins
        r_xy[i] = round( - L / 2 + (i - 0.5) * Delta_xy, digits = 6);
        println(Density_Distribution_File, "$(r_xy[i])\t$(round(g_x_Mean[i], digits = 6))\t$(round(g_x_Std[i], digits = 6))\t$(round(g_y_Mean[i], digits = 6))\t$(round(g_y_Std[i], digits = 6))")
    end
    @inbounds for i = 1:Nz_Bins
        r_z[i] = round( - h / 2 + (i - 0.5) * Delta_z, digits = 6);
        println(Density_Distribution_Z_File, "$(r_z[i])\t$(round(g_z_Mean[i], digits = 6))\t$(round(g_z_Std[i], digits = 6))")
    end
    close(Density_Distribution_File)
    close(Density_Distribution_Z_File)

    Distribution_X_Plot = Distributions_Plots("Density_Distribution_X", r_xy, g_x_Mean, g_x_Std, g_z_Max, h, Output_Route, true)
    Distribution_Y_Plot = Distributions_Plots("Density_Distribution_Y", r_xy, g_y_Mean, g_y_Std, g_z_Max, h, Output_Route, true)
    Distribution_Z_Plot = Distributions_Plots("Density_Distribution_Z", r_z, g_z_Mean, g_z_Std, g_z_Max, h, Output_Route, true)
    Distribution_X_Plot, Distribution_Y_Plot, Distribution_Z_Plot = Distribution_Plot_Unified("Density_Distribution", "Density", Distribution_X_Plot, Distribution_Y_Plot, Distribution_Z_Plot, Output_Route)

    ####################################################    CAVITY PROBABILITY AND MOVEMENT/INSERTION/DELETION OF MOLECULES  ################################################################
    Pc_Array = zeros(Float64, length(Pc) - 1)
    Pc_Density = zeros(Float64, length(Pc_Array))
    @inbounds for i in keys(Pc)
        i != 0 ? Pc_Array[i] = Pc[i] : nothing
    end
    Pc_File = open("$Output_Route/Cavity_Probability.dat", "w");
    println(Pc_File, "Density\tCavityProbability")
    @inbounds for i = 1:length(Pc_Array)
        Pc_Density[i] = i / V;
        println(Pc_File, "$(Pc_Density[i])\t$(round(Pc_Array[i], digits = 6))")
    end
    close(Pc_File)
    Cavity_Probability_Plot = plot(Pc_Density, Pc_Array, xlabel = "Density", ylabel = "Cavity Probability", ylim = (0, 1), legend = false, foreground_color_legend = false, background_color_legend = false, framestyle = :box, width = 3, legendfontsize = 20, guidefontsize = 20, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    savefig(Cavity_Probability_Plot, "$Output_Route/Cavity_Probability_Plot")

    ####################################################    ENERGY, AVERAGE ENERGY AND ENERGY HISTOGRAM  ################################################################
    Energy_Plot = plot(Energy_Array, xlabel = "N. Measurements", ylabel = L"U^* / N", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!(Energy_Plot, [Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    
    Energy_Histogram =  histogram(Energy_Array, normalize = true, xlabel = L"U^* / N", ylabel = "Normalized Frequency", legend = false, framestyle = :box, bins = 20, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    vline!([Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(Mean_Energy_Array[N_Measurements], STD_Energy_Array[N_Measurements]), width = 3, linecolor = :black)
    
    Mean_Energy_Plot = plot(Mean_Energy_Array, ribbon = STD_Energy_Array, xlabel = "N. Measurements", ylabel = L"\langle U^* / N \rangle", fillalpha = 0.2, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([Mean_Energy_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    
    Energy_Plots = plot(Energy_Plot, Energy_Histogram, Mean_Energy_Plot, layout = (@layout [a{0.3h} ; b c]))
    savefig(Energy_Plots, "$Output_Route/Energy_Plots")
    ####################################################    DENSITY, AVERAGE DENSITY AND DENSITY HISTOGRAM  ################################################################
    Density_Plot = plot(Density_Array, xlabel = "N. Measurements", ylabel = L"\rho^*", titlefontsize = 25, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!(Density_Plot, [Mean_Density_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    
    Density_Histogram =  histogram(Density_Array, normalize = true, xlabel = L"\rho^*", ylabel = " Normalized Frequency", legend = false, framestyle = :box, bins = 20, width = 3, guidefontsize = 22, tickfontsize = 18,  left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    vline!([Mean_Density_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    plot!(Normal(Mean_Density_Array[N_Measurements], STD_Density_Array[N_Measurements]), width = 3, linecolor = :black)
    
    Mean_Density_Plot = plot(Mean_Density_Array, ribbon = STD_Density_Array, xlabel = "N. Measurements", ylabel = L"\langle \rho^* \rangle", fillalpha = 0.2, legend = false, framestyle = :box, width = 3, guidefontsize = 22, tickfontsize = 18, left_margin = 5mm, bottom_margin = 5mm, top_margin = 5mm, widen = true, size = [1920, 1080], dpi = 300)
    hline!([Mean_Density_Array[N_Measurements]], color = :black, width = 2, linestyle = :dash)
    
    Density_Plots = plot(Density_Plot, Density_Histogram, Mean_Density_Plot, layout = (@layout [a{0.3h} ; b c]))
    savefig(Density_Plots, "$Output_Route/Density_Plots")
    ####################################################    POVRAY .pov AND .ini FILES GENERATION  ################################################################
    if PovRay
        Povray_ini(h, L, μ, T, N_Image - 1)
        Povray_Pov(h, L, μ, T, σ_p, σ_w)
        run(`povray $Output_Route/Positions/Pore_X_Axis.ini`)
        run(`povray $Output_Route/Positions/Pore_Z_Axis.ini`)
        run(`povray $Output_Route/Positions/Slab.pov +W1280 +H720`)
    end
end

function Movement(h::Type, L::Type, Beta::Type, Z_Displacement::Float64, Displacement::Float64, Energy::Float64, N_Displacement_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, σ_p::Type = 1., λ_p::Type = 1.5, σ_w::Type = 1., λ_w::Type = 1.5, R_Cut::Type = 3.) where {Type <: Real}
    i = rand(1:length(x));
    Energy_Old = Energy_Calculation(h, L, x[i], y[i], z[i], x, y, z, σ_p, λ_p, σ_w, λ_w);
    x_Old, y_Old, z_Old = x[i], y[i], z[i];
    x[i] += Displacement * (rand() - 0.5);
    y[i] += Displacement * (rand() - 0.5);
    z[i] += Z_Displacement * (rand() - 0.5);
    x[i] = PeriodicBoundaryConditions!(L, x[i]);
    y[i] = PeriodicBoundaryConditions!(L, y[i]);
    z[i] = PeriodicBoundaryConditions!(L, z[i]);
    Energy_New = Energy_Calculation(h, L, x[i], y[i], z[i], x, y, z, σ_p, λ_p, σ_w, λ_w);
    Delta_E = Energy_New - Energy_Old;
    if rand() < exp(-Beta * Delta_E)
        N_Displacement_Accepted += 1;
        Energy += Delta_E;
    else
        x[i], y[i], z[i] = x_Old, y_Old, z_Old;
    end
    return Energy, N_Displacement_Accepted
end

function Insertion(h::Type, L::Type, μ::Type, Beta::Type, Energy::Float64, N_Insertion_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, σ_p::Type = 1., λ_p::Type = 1.5, σ_w::Type = 1., λ_w::Type = 1.5, R_Cut::Type = 3.) where {Type <: Real}
    x_Insertion, y_Insertion, z_Insertion = L * (rand() - 0.5), L * (rand() - 0.5), (h + σ_w) * (rand() - 0.5)
    Energy_Insertion = Energy_Calculation(h, L, x_Insertion, y_Insertion, z_Insertion, x, y, z, σ_p, λ_p, σ_w, λ_w)
    if rand() < exp( Beta * (μ - Energy_Insertion) + log((L^2. * (h + σ_w)) / (length(x) + 1)) )
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion), append!(y, y_Insertion), append!(z, z_Insertion)
        Energy += Energy_Insertion;
    end
    return Energy, N_Insertion_Accepted
end

function Insertion_Mezei(h::Float64, L::Type, μ::Type, Beta::Type, Energy::Float64, N_Insertion_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, x_Insertion::Array{Float64, 1}, y_Insertion::Array{Float64, 1}, z_Insertion::Array{Float64, 1}, Pc::Dict{Int64, Float64}, σ_p::Type = 1., λ_p::Type = 1.5, σ_w::Type = 1., λ_w::Type = 1.5, R_Cut::Type = 3.) where {Type <: Real}
    i = rand(1:length(x_Insertion))
    Energy_Insertion = Energy_Calculation(h, L, x_Insertion[i], y_Insertion[i], z_Insertion[i], x, y, z, σ_p, λ_p, σ_w, λ_w);
    if rand() < ((L^2. * (h + σ_w)) * Pc[length(x)] / (length(x) + 1) ) * exp(Beta * (μ - Energy_Insertion))
        N_Insertion_Accepted += 1;
        append!(x, x_Insertion[i]), append!(y, y_Insertion[i]), append!(z, z_Insertion[i])
        Energy += Energy_Insertion;
    end
    return Energy, N_Insertion_Accepted
end

function Removal(h::Type, L::Type, μ::Type, Beta::Type, Energy::Float64, N_Removal_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, σ_p::Type = 1., λ_p::Type = 1.5, σ_w::Type = 1., λ_w::Type = 1.5, R_Cut::Type = 3.) where {Type <: Real}
    i = rand(1:length(x));
    Energy_Removal = Energy_Calculation(h, L, x[i], y[i], z[i], x, y, z, σ_p, λ_p, σ_w, λ_w);
    if rand() < exp( Beta * (Energy_Removal - μ) + log(length(x) / (L^2. * (h + σ_w))) )
        N_Removal_Accepted += 1;
        deleteat!(x, i), deleteat!(y, i), deleteat!(z, i)
        Energy -= Energy_Removal;
    end
    return Energy, N_Removal_Accepted
end

function Removal_Mezei(h::Float64, L::Type, μ::Type, Beta::Type, Energy::Float64, N_Removal_Accepted::Int64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, Pc_Interpolation::Float64, σ_p::Type = 1., λ_p::Type = 1.5, σ_w::Type = 1., λ_w::Type = 1.5, R_Cut::Type = 3.) where {Type <: Real}
    i = rand(1:length(x))
    Energy_Removal = Energy_Calculation(h, L, x[i], y[i], z[i], x, y, z, σ_p, λ_p, σ_w, λ_w);
    if rand() < ( length(x) / ((L^2. * (h + σ_w)) * Pc_Interpolation) ) * exp(Beta * (Energy_Removal - μ))
        N_Removal_Accepted += 1;
        deleteat!(x, i), deleteat!(y, i), deleteat!(z, i)
        Energy -= Energy_Removal;
    end
    return Energy, N_Removal_Accepted
end

function Energy_Calculation(h::Type, L::Type, rx::Float64, ry::Float64, rz::Float64, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, σ_p::Type = 1., λ_p::Type = 1.5, σ_w::Type = 1., λ_w::Type = 1.5, R_Cut::Type = 3.) where {Type <: Real}
    Energy = 0;
    ####################################################    ENERGY CONTRIBUTION FROM THE SLAB  ################################################################
    if abs(rz) > h / 2 - λ_w
        x_Index = ceil(rx + L / 2) / σ_w + 2;
        y_Index = ceil((ry + L / 2) / (2 * √3 * σ_w/2) + √3);
        for i_y = y_Index - 2:y_Index + 2, i_x = x_Index - 2:x_Index + 2

            x_Position = - L / 2 - 1.5σ_w + (i_x - 1) * σ_w;
            y_Position = - L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w;
            if x_Position < rx + λ_w && x_Position > rx - λ_w
                    if y_Position < ry + λ_w && y_Position > ry - λ_w
                        Delta_x, Delta_y, Delta_z = rx - x_Position, ry - y_Position, abs(rz) - (h / 2);
                        r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
                        Energy += U(r, σ_w, λ_w, 1.0);
                        Energy == Inf ? (return Energy) : nothing
                    end
            end

            x_Position = - L / 2 - 1.5σ_w + (i_x - 1) * σ_w + σ_w/2;
            y_Position = - L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w + √3/2 * σ_w;
            if x_Position < rx + λ_w && x_Position > rx - λ_w
                    if y_Position < ry + λ_w && y_Position > ry - λ_w
                        Delta_x, Delta_y, Delta_z = rx - x_Position, ry - y_Position, abs(rz) - (h / 2);
                        r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
                        Energy += U(r, σ_w, λ_w, 1.0);
                        Energy == Inf ? (return Energy) : nothing
                    end
            end
        end
    end
    ####################################################    ENERGY CONTRIBUTION BETWEEN PARTICLES  ################################################################
    @inbounds for i = 1:length(x)
        Delta_x, Delta_y, Delta_z = rx - x[i], ry - y[i], rz - z[i];
        Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
        Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
        r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
        if r != 0. && r < R_Cut
            Energy += U(r, σ_p, λ_p, 1.0);
        end
        Energy == Inf ? (return Energy) : nothing
    end

    return Energy
end

function U(r::Float64, σ::Type = 1., λ::Type = 1.5, e::Type = 1.) where {Type <: Real}
    r <= σ ? (return Inf) : r <= λ ? (return -e) : (return 0)
end

function PeriodicBoundaryConditions!(L::Type, x::Float64) where {Type <: Real}
    return x - L * round(x / L)
end

function Random_Excluded_Volume(Equilibrium::Bool, h::Type, L::Type, Pc::Dict{Int64, Float64}, Pc_Sum::Dict{Int64, Float64}, Pc_N::Dict{Int64, Int64}, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, σ_p::Type = 1., σ_w::Type = 1.) where {Type <: Real}
    N_in = 0;
    Equilibrium ? N_Random = 200 : N_Random = 1000
    x_Insertion, y_Insertion, z_Insertion = Float64[], Float64[], Float64[];
    if !haskey(Pc, length(x))
            Equilibrium = false;
            N_Random = 1000;
    end
    @inbounds for i = 1:N_Random
        Control = false;
        x_V, y_V, z_V = L * (rand() - 0.5), L * (rand() - 0.5), h * (rand() - 0.5);
        for j = 1:length(x)
            Delta_x, Delta_y, Delta_z = x_V - x[j], y_V - y[j], z_V - z[j];
            Delta_x = PeriodicBoundaryConditions!(L, Delta_x);
            Delta_y = PeriodicBoundaryConditions!(L, Delta_y);
            r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
            if r < σ_p
                Control = true;
                break
            end
        end
        if abs(z_V) >= (h - σ_w) / 2
            if Control == false
                for i_y = 1:ceil(L / (√3 * σ_w) + √3), i_x = 1:ceil(L / σ_w + 3 + 1)
                    Delta_x = x_V - (- L / 2 - 1.5σ_w + (i_x - 1) * σ_w);
                    Delta_y = y_V - (- L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w);
                    Delta_z = abs(z_V) - ((h + σ_w) / 2)
                    r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
                    if r < σ_w
                        Control = true;
                        break
                    end
                    Delta_x = x_V - (- L / 2 - 1.5σ_w + (i_x - 1) * σ_w + σ_w/2);
                    Delta_y = y_V - (- L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w + √3 * σ_w/2);
                    Delta_z = abs(z_V) - ((h + σ_w) / 2)
                    r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
                    if r < σ_w
                        Control = true;
                        break
                    end
                end
            end
        end
        if Control == false
            N_in += 1;
            append!(x_Insertion, x_V), append!(y_Insertion, y_V), append!(z_Insertion, z_V)
        end
    end
    if Equilibrium == false
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
    end
    return Pc, Pc_Sum, Pc_N, x_Insertion, y_Insertion, z_Insertion
end

function Grid_Excluded_Volume(h::Type, L::Type, Pc_Grid::Dict{Int64, Float64}, Pc_Grid_Sum::Dict{Int64, Float64}, Pc_Grid_N::Dict{Int64, Int64}, x::Array{Float64, 1}, y::Array{Float64, 1}, z::Array{Float64, 1}, σ_p::Type = 1., σ_w::Type = 1.) where {Type <: Real}
    Grid, Sep_xy, Sep_z, EVMPS = 10, L / (Grid + 1), h / (Grid + 1), zeros(Grid, Grid, Grid);
    @inbounds for i_x = 1:Grid
        @inbounds for i_y = 1:Grid
            @inbounds for i_z = 1:Grid
                Control = false;
                x_Grid, y_Grid, z_Grid = i_x * Sep_xy - L / 2., i_y * Sep_xy - L / 2., i_z * Sep_z - (h + σ_w) / 2.;
                @inbounds for i = 1:length(x)
                    Delta_x, Delta_y, Delta_z = x_Grid - x[i], y_Grid - y[i], z_Grid - z[i];
                    r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
                    if r < σ_p
                        Control, EVMPS[i_x, i_y, i_z] = true, 1;
                        break
                    end
                end
                if Control == false
                    for i_y = 1:ceil(L / (√3 * σ_w) + √3), i_x = 1:ceil(L / σ_w + 3 + 1)
                        Delta_x = x_Grid - (- L / 2 - 1.5σ_w + (i_x - 1) * σ_w);
                        Delta_y = y_Grid - (- L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w);
                        Delta_z = abs(z_Grid) - ((h + σ_w) / 2);
                        r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
                        if r < σ_w
                            Control, EVMPS[i_x, i_y, i_z] = true, 1;
                            break
                        end
                        Delta_x = x_Grid - (- L / 2 - 1.5σ_w + (i_x - 1) * σ_w + σ_w/2);
                        Delta_y = y_Grid - (- L / 2 - 1.5σ_w + (i_y - 1) * √3 * σ_w + √3 * σ_w/2);
                        Delta_z = abs(z_Grid) - ((h + σ_w) / 2)
                        r = sqrt(Delta_x^2 + Delta_y^2 + Delta_z^2);
                        if r < σ_w
                            Control, EVMPS[i_x, i_y, i_z] = true, 1;
                            break
                        end
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

function Distribution(N_Bins::Int64, L::Type, x::Array{Float64, 1}) where {Type <: Real}
    Delta = L / N_Bins;
    g = zeros(Float64, N_Bins);
    for i = 1:length(x)
        l = convert(Int64, ceil( (L / 2. - x[i]) / Delta) )
        g[l] += 1
    end
    return g
end

function Povray_Pov(h::Type, L::Type, μ::Type, T::Type, σ_p::Type = 1., σ_w::Type = 1.) where {Type <: Real}
    Particle_r, Particle_g, Particle_b, Particle_t = 174/255, 214/255, 241/255, 0;
    Slab_Attractive_r, Slab_Attractive_g, Slab_Attractive_b, Slab_Attractive_t = 241/255, 196/255, 15/255, 0; 
    Slab_Neutral_r, Slab_Neutral_g, Slab_Neutral_b, Slab_Neutral_t  = 142/255, 68/255, 173/255, 0; 
    Slab_Repulsive_r, Slab_Repulsive_g, Slab_Repulsive_b, Slab_Repulsive_t = 183/255, 28/255, 28/255, 0;
    Wall_r, Wall_g, Wall_b, Wall_t = 174/255, 214/255, 241/255, 0.5;
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(μ, digits = 3))/h_$(h)/Positions"
    ###################################################      X AXIS VIEW     ################################################################
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
    println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w)\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\nsphere {\n\t\t<rx, ry, rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n#end")
    println(Pov_File, "#macro Wall(L, h)\n\tunion{\n\t\ttriangle {\n\t\t\t<-L / 2, -L / 2, h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t\ttriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t}
        union {\n\t\ttriangle {\n\t\t\t<L / 2, -L / 2, h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n\n\t
        union {\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <-L / 2, L / 2, -h /2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\ntriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <L / 2, L / 2, h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n#end")
    println(Pov_File, "#declare L = $L;\n#declare h = $h;")
    println(Pov_File, """#fopen File_Positions concat("$Output_Route/Pos_", str(clock, 1, 0), ".xyz") read""")
    println(Pov_File, "\t#while (defined( File_Positions ))\n\t\t#read (File_Positions, rx, ry, rz)\n\t\tParticle(rx, ry, -rz)\n\t\t#declare PBC = false;\n\t\t#if (rx > (L - 1) / 2)\n\t\t\t#declare rx = rx - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (rx < -(L - 1) / 2)\n\t\t\t#declare rx = rx + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t
        #if (ry > (L - 1) / 2)\n\t\t\t#declare ry = ry - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (ry < -(L - 1) / 2)\n\t\t\t#declare ry = ry + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (PBC)\n\t\t\tParticle(rx, ry, -rz)\n\t\t#end\n\t#end\n#fclose File_Positions")
    println(Pov_File, """#fopen File_Slab "$Output_Route/Slab.xyz" read""")
    println(Pov_File, "#while (defined (File_Slab)) \n\t#read (File_Slab, rx, ry, rz)\n\tSlab(rx, ry, rz, $(σ_w/2))\n#end\n#fclose File_Slab")
    println(Pov_File, "Wall(L, h)")
    close(Pov_File)
    ###################################################      Z AXIS VIEW     ################################################################
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
    println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w)\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n#end")
    println(Pov_File, "#macro Wall(L, h)\n\tunion{\n\t\ttriangle {\n\t\t\t<-L / 2, -L / 2, h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t\ttriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <-L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t}
        union {\n\t\ttriangle {\n\t\t\t<L / 2, -L / 2, h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <L / 2, -L / 2, -h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n\n\t
        union {\n\t\ttriangle {\n\t\t\t<L / 2, L / 2, -h /2>, <L / 2, L / 2, h / 2>, <-L / 2, L / 2, -h /2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\ntriangle {\n\t\t\t<-L / 2, L / 2, -h /2>, <-L / 2, L / 2, h / 2>, <L / 2, L / 2, h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n
        union {\n\t\ttriangle {\n\t\t\t<L / 2, -L / 2, -h /2>, <L / 2, -L / 2, h / 2>, <-L / 2, -L / 2, -h /2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\ntriangle {\n\t\t\t<-L / 2, -L / 2, -h /2>, <-L / 2, -L / 2, h / 2>, <L / 2, -L / 2, h / 2>\n\t\t\tpigment { rgbt <$Wall_r, $Wall_g, $Wall_b, $Wall_t> }\n\t\t}\n\tno_shadow}\n#end")
    println(Pov_File, "#declare L = $L;\n#declare h = $h;")
    println(Pov_File, """#fopen File_Positions concat("$Output_Route/Pos_", str(clock, 1, 0), ".xyz") read""")
    println(Pov_File, "\t#while (defined( File_Positions ))\n\t\t#read (File_Positions, rx, ry, rz)\n\t\tParticle(rx, ry, -rz)\n\t\t#declare PBC = false;\n\t\t#if (rx > (L - 1) / 2)\n\t\t\t#declare rx = rx - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (rx < -(L - 1) / 2)\n\t\t\t#declare rx = rx + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t
        #if (ry > (L - 1) / 2)\n\t\t\t#declare ry = ry - L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (ry < -(L - 1) / 2)\n\t\t\t#declare ry = ry + L;\n\t\t\t#declare PBC = true;\n\t\t#end\n\t\t#if (PBC)\n\t\t\tParticle(rx, ry, -rz)\n\t\t#end\n\t#end\n#fclose File_Positions")
    println(Pov_File, """#fopen File_Slab "$Output_Route/Slab.xyz" read""")
    println(Pov_File, "#while (defined (File_Slab)) \n\t#read (File_Slab, rx, ry, rz)\n\tSlab(rx, ry, rz, $(σ_w/2))\n#end\n#fclose File_Slab")
    println(Pov_File, "Wall(L, h)")
    close(Pov_File)
    ##################################################      WALL DECORATION     ################################################################
    Pov_File = open("$Output_Route/Slab.pov", "w");
    println(Pov_File, "global_settings {\n\tambient_light rgb <0.2, 0.2, 0.2>\tmax_trace_level 15\n}\n")
    println(Pov_File, "background { color rgb <1, 1, 1> }\n")
    println(Pov_File, "#default { finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.7 phong 1} }\n")
    println(Pov_File, "camera {\n\tperspective\n\tlocation <0, 0, $(1.3L)>\n\tlook_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(-5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(+5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(-5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(+5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w)\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n#end")
    println(Pov_File, """#fopen File_Slab "$Output_Route/Slab.xyz" read""")
    println(Pov_File, "#while (defined (File_Slab)) \n\t#read (File_Slab, rx, ry, rz)\n\tSlab(rx, ry, 0, $(σ_w/2))\n#end\n#fclose File_Slab")
    close(Pov_File)
end

function Povray_ini(h::Type, L::Type, μ::Type, T::Type, Frames::Int64) where {Type <: Real}
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(μ, digits = 3))/h_$(h)/Positions"
    Ini_File = open("$Output_Route/Pore_X_Axis.ini", "w");
    println(Ini_File, "Input_File_Name = $Output_Route/Pore_X_Axis.pov\nOutput_File_Name = $Output_Route/\n+W1280 +H720\nInitial_Frame = 1\nFinal_Frame = $Frames\nInitial_Clock = 1\nFinal_Clock = $Frames\nCyclic_Animation = off")
    close(Ini_File)

    Ini_File = open("$Output_Route/Pore_Z_Axis.ini", "w");
    println(Ini_File, "Input_File_Name = $Output_Route/Pore_Z_Axis.pov\nOutput_File_Name = $Output_Route/\n+W1280 +H720\nInitial_Frame = 1\nFinal_Frame = $Frames\nInitial_Clock = 1\nFinal_Clock = $Frames\nCyclic_Animation = off")
    close(Ini_File)
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

Args = Parse_Commandline()
println("Parsed args:")
for (arg,val) in Args
    println("  $arg  =>  $val")
end
@time Confined_GrandCanonical_MonteCarlo(Args["μ"], Args["T"], Args["h"], Args["ρ_Bulk"], Args["Configurations"], Args["PovRay"])