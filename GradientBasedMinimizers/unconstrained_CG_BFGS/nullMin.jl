# Julia implementation of the nullspace formulation for minimizing solution phase models

using StaticArrays;
using LinearAlgebra;
using BenchmarkTools
using Statistics
using Plots

# change path to the location of this file
cd(@__DIR__)

include("functions/sys_info.jl")
include("functions/figures/figures.jl")

"""
    select_ph!(ph::Symbol)

Selects and includes the appropriate phase and gradient method files based on the provided phase symbol.
"""
function select_ph!(ph::Symbol)
    if ph == :hb
        include("functions/phases/hb.jl")
        include("functions/gradient_method.jl")
    elseif ph == :cpx
        include("functions/phases/cpx.jl")
        include("functions/gradient_method.jl")
    elseif ph == :spn
        include("functions/phases/spn.jl")
        include("functions/gradient_method.jl")
    else
        error("Unknown phase: $ph")
    end
end

#choose the reference Gibbs hyperplane: em, pc, eq (endmember only, endmember+pseudocompounds, equilibrium):
gam  = "eq";
#choose gradient based minimization method: cg, BFGS (conjugate-gradient, BFGS)
gopt = :BFGS;
gv = init_sys_infos();

# select phase (:hb, :cpx, :spn)
ph = :spn;
select_ph!(ph);

#choose test (only for spinel, choose 1 or 2. 2 -> solvus)
test = 1;

#Initialize system informations
ph = init_phase(gam,test);              #initialize phase to minimize
gm = init_gradient_methods(ph);    #initialize gradient method

if gopt == :BFGS
    null_min! = ((gv,ph,gm)) -> null_min_BFGS!(gv,ph,gm);
elseif gopt == :CG
    null_min! = ((gv,ph,gm)) -> null_min_CG!(gv,ph,gm);
end

# retrieve the nullspace
update_Nullspace!(gm,ph);

# generate figures
generate_min_time_vs_distance_figure(ph, gv, gm)
generate_min_time_vs_normDeltaGamma_figure(ph, gv, gm)
create_box_plot()

if test == 2
    generate_figure_solvus(ph, gv, gm)
end