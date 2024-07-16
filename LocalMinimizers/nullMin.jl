# Julia implementation of the nullspace formulation for minimizing solution phase models

using StaticArrays;
using LinearAlgebra;
using BenchmarkTools
using Statistics
using Test
using Plots

include("sys_info.jl")

#=
    Choose the solution phase to optimize here
=#
include("hb.jl")
# include("cpx.jl")            
# include("spn.jl")


include("gradient_method.jl")           #import gradient based minimization methods

gam  = "eq";                            #choose the reference Gibbs hyperplane: em, pc, eq (endmember only, endmember+pseudocompounds, equilibrium):  
gopt = "BFGS";                          #choose gradient based minimization method: cg, BFGS (conjugate-gradient, BFGS)
test = 1                                #choose test (only for spinel, choose 1 or 2. 2 -> solvus)

gv = init_sys_infos();                  #Initialize system informations
ph = init_phase(gam,test);              #initialize phase to minimize
gm = init_gradient_methods(gopt,ph);    #initialize gradient method


if gopt == "BFGS"
    null_min! = ((gv,ph,gm)) -> null_min_BFGS!(gv,ph,gm);
elseif gopt == "CG"
    null_min! = ((gv,ph,gm)) -> null_min_CG!(gv,ph,gm);
end

# retrieve the nullspace
update_Nullspace!(gm,ph);

# run test point
ph.ig        .= get_ig(ph.pc_list[1,:]);

@benchmark null_min!(gv,ph,gm); 

# Generate figures
if 1==1
    include("Fig_minTime_vs_distance.jl")
    include("Fig_minTime_vs_normDeltaGamma.jl")
    # include("Fig_solvus.jl")                            #select test 2 above, test 1 is without solvus for spinel
    include("Fig_mean_min_time.jl")
end


