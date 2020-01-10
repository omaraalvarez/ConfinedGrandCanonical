using Plots;
using CSV;

Percentage_Array = [0, 10, 20, 30, 40, 50, 60, 70, 80, 100];
T = 1.5;
ChemPot = -3.0;
p = plot(foreground_color_legend = false, legend = :bottomright, xlabel = "Slit Separation", ylabel= "< Density > [Unitless]", size = [1920, 1080])
Route = pwd() * "/Output/T_$(T)/ChemPot_$ChemPot"
for Percentage in Percentage_Array
    a = CSV.read("$Route/Patch_$Percentage%/Density_T_$(T).dat", delim = "\t");
    plot!(a.h, a.Density, yerror = a.ErrorDensity, tickfontsize = 12, guidefontsize = 12, legendfontsize = 12, width = 3, label = "$Percentage% Attractive Area")
end
savefig(p, "$Route/Adsorption_Plot")