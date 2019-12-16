function Povray_Pov(h::Float64, L::Float64, ChemPot::Float64, T::Float64, σ_w::Float64, λ_w::Float64)
    Particle_r, Particle_g, Particle_b, Particle_t = 174/255, 214/255, 241/255, 0;
    Slab_Attractive_r, Slab_Attractive_g, Slab_Attractive_b, Slab_Attractive_t = 26/255, 35/255, 126/255, 0; 
    Slab_Neutral_r, Slab_Neutral_g, Slab_Neutral_b, Slab_Neutral_t  = 77/255, 86/255, 86/255, 0; 
    Slab_Repulsive_r, Slab_Repulsive_g, Slab_Repulsive_b, Slab_Repulsive_t = 183/255, 28/255, 28/255, 0;
    Wall_r, Wall_g, Wall_b, Wall_t = 174/255, 214/255, 241/255, 0.5;

    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/lambdaw_$(λ_w)/h_$(h)/Positions"
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
    if λ_w > 1.
        println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n\t\n\t\n\tsphere {\n\t\t<rx, ry, rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n#end")
    else
        println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Neutral_r, $Slab_Neutral_g, $Slab_Neutral_b, $Slab_Neutral_t> }\n\tno_shadow}\n\tsphere {\n\t\t<rx, ry, rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Neutral_r, $Slab_Neutral_g, $Slab_Neutral_b, $Slab_Neutral_t> }\n\tno_shadow}\n#end")
    end    
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

function Povray_Pov_Z_Axis(h::Float64, L::Float64, ChemPot::Float64, T::Float64, σ_w::Float64, λ_w::Float64)
    Particle_r, Particle_g, Particle_b, Particle_t = 174/255, 214/255, 241/255, 0;
    Slab_Attractive_r, Slab_Attractive_g, Slab_Attractive_b, Slab_Attractive_t = 26/255, 35/255, 126/255, 0; 
    Slab_Neutral_r, Slab_Neutral_g, Slab_Neutral_b, Slab_Neutral_t  = 77/255, 86/255, 86/255, 0; 
    Slab_Repulsive_r, Slab_Repulsive_g, Slab_Repulsive_b, Slab_Repulsive_t = 183/255, 28/255, 28/255, 0;
    Wall_r, Wall_g, Wall_b, Wall_t = 174/255, 214/255, 241/255, 0.5;
    
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/lambdaw_$(λ_w)/h_$(h)/Positions"
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
    if λ_w > 1.
        println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n\t\n#end")
    else
        println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Neutral_r, $Slab_Neutral_g, $Slab_Neutral_b, $Slab_Neutral_t> }\n\tno_shadow}\n\t\n#end")
    end    
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

function Povray_Slab(L::Float64, ChemPot::Float64, T::Float64, σ_w::Float64, λ_w::Float64)
    Slab_Attractive_r, Slab_Attractive_g, Slab_Attractive_b, Slab_Attractive_t = 26/255, 35/255, 126/255, 0; 
    Slab_Neutral_r, Slab_Neutral_g, Slab_Neutral_b, Slab_Neutral_t  = 77/255, 86/255, 86/255, 0; 
    Slab_Repulsive_r, Slab_Repulsive_g, Slab_Repulsive_b, Slab_Repulsive_t = 183/255, 28/255, 28/255, 0;
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/lambdaw_$(λ_w)"
    Pov_File = open("$Output_Route/Pore_Slab.pov", "w");
    println(Pov_File, "global_settings {\n\tambient_light rgb <0.2, 0.2, 0.2>\tmax_trace_level 15\n}\n")
    println(Pov_File, "background { color rgb <1, 1, 1> }\n")
    println(Pov_File, "#default { finish {ambient .8 diffuse 1 specular 1 roughness .005 metallic 0.7 phong 1} }\n")
    println(Pov_File, "camera {\n\tperspective\n\tlocation <0, 0, $(1.3L)>\n\tlook_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(-5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<0, $(+5L), 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(-5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    println(Pov_File, "light_source {\n\t<$(+5L), 0, 0>\n\tcolor rgb <0.3, 0.3, 0.3>\n\tfade_distance $(10L)\n\tfade_power 0\n\tparallel\n\tpoint_at <0, 0, 0>\n}\n")
    if λ_w > 1.
        println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Attractive_r, $Slab_Attractive_g, $Slab_Attractive_b, $Slab_Attractive_t> }\n\tno_shadow}\n\t\n#end")
    else
        println(Pov_File, "#macro Slab(rx, ry, rz, sigma_w, lambda_w)\n\t\n\tsphere {\n\t\t<rx, ry, -rz>, sigma_w\n\t\tpigment {rgbt <$Slab_Neutral_r, $Slab_Neutral_g, $Slab_Neutral_b, $Slab_Neutral_t> }\n\tno_shadow}\n\t\n#end")
    end
    println(Pov_File, """#fopen File_Slab "$Output_Route/Slab.xyz" read""")
    println(Pov_File, "#while (defined (File_Slab)) \n\t#read (File_Slab, rx, ry, rz)\n\tSlab(rx, ry, 0, $σ_w, $λ_w)\n#end\n#fclose File_Slab")
    close(Pov_File)
end

function Povray_ini(h::Float64, ChemPot::Float64, T::Float64, λ_w::Float64, Frames::Int64)
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/lambdaw_$(λ_w)/h_$(h)/Positions"
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

function Povray_ini_Z_Axis(h::Float64, ChemPot::Float64, T::Float64, λ_w::Float64, Frames::Int64)
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/lambdaw_$(λ_w)/h_$(h)/Positions"
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

function Povray_Images()
    H = range(2., 10., step = 1.);
    T = 1.25
    ChemPot = -2.
    L = 20.
    λ_w = 1.5
    Output_Route = pwd() * "/Output/T_$(round(T, digits = 2))/ChemPot_$(round(ChemPot, digits = 2))/lambdaw_$(λ_w)"

    for h in H
        Povray_Pov(h, L, ChemPot, T, 0.5, λ_w)
        Povray_ini(h, ChemPot, T, λ_w, 20)
        run(`povray $Output_Route/h_$(h)/Positions/Pore_Animation.ini`)
        Povray_Pov_Z_Axis(h, L, ChemPot, T, 0.5, λ_w)
        Povray_ini_Z_Axis(h, ChemPot, T, λ_w, 20)
        run(`povray $Output_Route/h_$(h)/Positions/Pore_Z_Axis_Animation.ini`)
    end
    Povray_Slab(L, ChemPot, T, 0.5, λ_w)
    run(`povray $Output_Route/Pore_Slab.pov `)
end

Povray_Images()