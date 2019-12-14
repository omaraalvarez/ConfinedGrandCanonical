using Plots;
using CSV;

T_Array = [0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0]
ChemPot = [-3., -2., -1., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.]

Average_Density_Plot = plot(background_color_legend = false, foreground_color_legend = false, legend = :bottomright, xlabel = "Chemical Potential", ylabel= "< Density > [Unitless]", width = 3, size = [1200, 800])
for T in T_Array
    Density_Route = pwd() * "/Output_Julia/T_$(round(T, digits = 2))"
    Density = CSV.read("$Density_Route/Density_T_$T.dat", delim = "\t")
    plot!(Density.ChemPot, Density.Density, yerror = Density.ErrorDensity, label = "T = $T", width = 3)
end
Plot_Route = pwd() * "/Output_Julia"
savefig(Average_Density_Plot, "$Plot_Route/Density")

#for T in T_Array
#    Density_Route = pwd() * "/Output_Julia/T_$(round(T, digits = 2))"
#    p = Array{Any, 1}(undef, length(ChemPot))
#    i = 1
#    for μ in ChemPot
#        Chem_Route = pwd() * "/Output_Julia/T_$(round(T, digits = 2))/ChemPot_$(round(μ, digits = 2))"
#        Radial = CSV.read("$Chem_Route/Radial_Distribution.dat", delim = "\t", header = false)
##        #plot!(Radial.Column1, Radial.Column2, label = "Chemical Potential = $μ", width = 3)
#        p[i] = plot((Radial.Column1, Radial.Column2), title = "Chemical Potential = $μ", grid = false, xlabel = "Distance [r]", ylabel= "g(r) [Unitless]", width = 3, legend = false)
#        hline!([1.], style = :dash, width = 3, color = :black)
#        i += 1
#    end
##    #hline!([1.], style = :dash, width = 3, color = :black, label = 1.0)
#    y = ones(3)
#    title = Plots.scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], Plots.text("T = $T", :black, 25)), axis=false, grid = false, leg=false,size=(1920,50))
#    Radial_Plot = plot(p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15], p[16],layout = (4, 4), size = [1920, 1080])
#    Final_Plot = plot(title, Radial_Plot, layout = grid(2, 1, heights = [0.025, 0.975]))
#    savefig(Final_Plot, "$Density_Route/Radial_Distribution_$T")
#end