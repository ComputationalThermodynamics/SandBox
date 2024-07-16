#=
    NR 5.07.24
    Toy code to test the performances of the C implementation of the local minimizer using NLopt
    By default the SLSQP algorithm is active in MAGEMin_C. 
    To try the CCSAQ algorithm, you need to download MAGEMin and the C code using CCSAQ instead of SLSQP
=#

include("functions/functions.jl")

using LinearAlgebra, Statistics, Plots
using MAGEMin_C, BenchmarkTools

# select phase (hb, cpx, spn)
ph = "spn"

if ph == "hb"
    select_ph!=  (()) -> select_hb!();
elseif ph == "cpx"
    select_ph! =  (()) -> select_cpx!();
elseif ph == "spn"
    select_ph! =  (()) -> select_spn!();
end

# Initialize MAGEMin
data        =   Initialize_MAGEMin("ig", verbose=true);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);

# Get the information for the minimization tests
gamma, pc, P, T, ph_id, ig, gamma_eq, sol = select_ph!()

# here we perform a minimization to initialize all needed variables
out         =   point_wise_minimization(P,T, data);

# test the method
Fig_minTime_vs_normDeltaGamma(ph, gamma, pc, P, T, ph_id, ig, gamma_eq, sol)
Fig_mean_min_time(ph, gamma, pc, P, T, ph_id, ig, gamma_eq, sol)


