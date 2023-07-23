<img src="./pics/Figure_LP_vs_PGE.png" alt="drawing" width="640" alt="centered image"/>

# Thermodynamic Sandbox

Here you can find script(s) dealing with thermodynamic-related problems

## Nullspace minimizer

> #### *NullspaceMinimization/nullMin.jl*
>
> - Julia scripts to use nullspace optimization approach to minimize individual site-fraction based solution phases.
> - The scripts include spinel, clinopyroxene and amphibole as formulated in Holland et al., 2018.
> - Gibbs hyperplane for each phase have been computed using MAGEMin and the igneous database.

	1. 12 kbar and 1100 °C for spinel
	2. 3.26 kbar and 906.25 °C for the spinel solvus test
	2. 12 kbar and 1100 °C for clinopyroxene
	3. 5 kbar and 650 °C for amphibole


## Gibbs free energy minimization

> #### *PhaseEquilibriumMinimizer/MAGEMin_PGE_and_LP.m*
>
> - MATLAB script used in Riel et al., (2022) [Geochemistry, Geophysics, Geosystems, 23, 7] to present a simplified application of MAGEMin minimization approach (including feldspar, quartz and sillimanite) and compare it with linear programming. 
> - The minimization is conducted in the NCKAS chemical system at 0.3 GPa and 873.15 K.
> - The Plagioclase model is taken after Holland et al., 2021
> - The Linear programming approach described in de Capitani & Brown (1987)

