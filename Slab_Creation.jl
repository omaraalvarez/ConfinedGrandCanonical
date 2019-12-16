Output_Route = pwd()
σ_w = 0.5
L = 20.
N_Slab = 2 * convert(Int64, ceil(L / (2σ_w) + 3 + 1) * ceil(L / (√3 * 2σ_w) + √3))
File_Slabs_Avogadro = open("$Output_Route/Slab_Avogadro.xyz", "w+")
println(File_Slabs_Avogadro, "$N_Slab\n")


################ PARCHE CIRCULAR ########################
#for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
#    r = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)^2 + (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)^2)
#    if r < 5.
#        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
#    else
#        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
#    end
#    r = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)^2 + (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)^2 )
#    if r < 5.
#        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
#    else
#        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
#    end
#end

##################### PARCHE CUADRADO #######################
#for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
#    if - L / 2 - 3σ_w + (i_x - 1) * 2σ_w > -L / 4 && - L / 2 - 3σ_w + (i_x - 1) * 2σ_w < L / 4 && - L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w > -L /4 && - L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w < L /4
#        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
#    else
#        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
#    end
#    if - L / 2 - 3σ_w + (i_x - 1) * 2σ_w > -L / 4 && - L / 2 - 3σ_w + (i_x - 1) * 2σ_w < L / 4 && - L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w > -L /4 && - L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w < L /4
#        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
#    else
#        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
#    end
#end

################### PARCHE ALTERNADO ###########################
#for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
#    if isodd(convert(Int64, i_x)) && isodd(convert(Int64, i_y))
#        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
#    else
#        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
#    end
#    if isodd(convert(Int64, i_x)) && isodd(convert(Int64, i_y))
#        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
#    else
#        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
#    end
#end

######################### PARCHE ALTERNADO 2 ####################
#for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
#    if iseven(convert(Int64, i_x)) && isodd(convert(Int64, i_y))
#        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
#    else
#        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
#    end
#    if isodd(convert(Int64, i_x)) && iseven(convert(Int64, i_y))
#        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
#    else
#        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
#    end
#end

for i_y = 1:ceil(L / (√3 * 2σ_w) + √3), i_x = 1:ceil(L / (2σ_w) + 3 + 1)
    r = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)^2 + (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)^2)
    if (r > 9 && r < 10) || (r > 7 && r < 8) || (r > 5 && r < 6) || (r > 3 && r < 4) || (r > 1 && r < 2)
        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
    else
        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)\t0.")
    end
    r = sqrt( (- L / 2 - 3σ_w + (i_x - 1) * 2σ_w)^2 + (- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w)^2)
    if (r > 9 && r < 10) || (r > 7 && r < 8) || (r > 5 && r < 6) || (r > 3 && r < 4) || (r > 1 && r < 2)        println(File_Slabs_Avogadro, "Li\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
    else
        println(File_Slabs_Avogadro, "H\t$(- L / 2 - 3σ_w + (i_x - 1) * 2σ_w + σ_w)\t$(- L / 2 - 3σ_w + (i_y - 1) * √3 * 2σ_w + √3 * σ_w)\t0.")
    end
end

close(File_Slabs_Avogadro)
