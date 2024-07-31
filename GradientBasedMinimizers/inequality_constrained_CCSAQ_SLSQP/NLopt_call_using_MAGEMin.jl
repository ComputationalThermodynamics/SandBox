#=
    NR 5.07.24
    Toy code to test the performances of the C implementation of the local minimizer using NLopt
    By default the SLSQP algorithm is active in MAGEMin_C. 
    To try the CCSAQ algorithm, you need to download MAGEMin and the C code using CCSAQ instead of SLSQP
=#

using LinearAlgebra, Statistics, Plots
using MAGEMin_C, BenchmarkTools

# change path to the location of this file
cd(@__DIR__)

include("functions/functions.jl")
include("functions/figures.jl")

# select phase (:hb, :cpx or :spn)
ph = :spn

# Initialize MAGEMin
data        =   Initialize_MAGEMin("ig", verbose=true);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);

# Get the information for the minimization tests
gamma, pc, P, T, ph_id, ig, gamma_eq, sol = select_ph(ph)

# here we perform a minimization to initialize all needed variables
out         =   point_wise_minimization(P,T, data);

# test the method and create plots
Fig_minTime_vs_normDeltaGamma(ph, gamma, pc, ph_id, ig, sol)
Fig_mean_min_time(ph, gamma, ph_id, ig, gamma_eq)
