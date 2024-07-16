#=
    NR 5.07.24
    Toy code to test the performances of the C implementation of the local minimizer using NLopt
=#
using LinearAlgebra;
using Statistics
using Plots

"""
    testjl(pc, gamma, ph_id, gv, DB)

    Test the Julia implementation of the local minimizer using NLopt

    # Arguments
    - pc :: Vector{Float64} : initial guess
    - gamma :: Vector{Float64} : initial chemical potentials
    - ph_id :: Int64 : phase id
    - gv :: LibMAGEMin.gv : system informations
    - DB :: LibMAGEMin.DB : database

    # Returns
    - xsol :: Vector{Float64} : optimized solution
"""
function testjl(pc      :: Vector{Float64},
                gamma   :: Vector{Float64},
                ph_id   :: Int64,
                gv, DB  )

    gv.maxeval   		= 1024;	
    gv.obj_tol   		= 1e-6;	

    # set the chemical potential Gamma

    # set the initial guess (xeos)
    SS_ref_db   = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

    unsafe_copyto!(gv.gam_tot, pointer(gamma), gv.len_ox)

    SS_ref_db[ph_id] = LibMAGEMin.rotate_hyperplane(gv,SS_ref_db[ph_id])
    SS_ref_db[ph_id] = LibMAGEMin.restrict_SS_HyperVolume(	gv,SS_ref_db[ph_id],1.0	);

    unsafe_copyto!(SS_ref_db[ph_id].iguess, pointer(pc), SS_ref_db[ph_id].n_xeos)

    if ph_id == 1
        SS_ref_db[ph_id] = LibMAGEMin.NLopt_opt_ig_spn_function(gv,SS_ref_db[ph_id])
    elseif ph_id == 4
        SS_ref_db[ph_id] = LibMAGEMin.NLopt_opt_ig_cpx_function(gv,SS_ref_db[ph_id])
    elseif ph_id == 7
        SS_ref_db[ph_id] = LibMAGEMin.NLopt_opt_ig_hb_function(gv,SS_ref_db[ph_id])
    end

    x   = unsafe_wrap(Vector{Float64},SS_ref_db[ph_id].xeos,SS_ref_db[ph_id].n_xeos)

    return xsol = copy(x)
end


"""
    testjl_rot(ig, gamma, ph_id, gv, DB)

    Test the Julia implementation of the local minimizer using NLopt

    # Arguments
    - ig :: Vector{Float64} : initial guess
    - gamma :: Vector{Float64} : initial chemical potentials
    - ph_id :: Int64 : phase id
    - gv :: LibMAGEMin.gv : system informations
    - DB :: LibMAGEMin.DB : database

    # Returns
    - xsol :: Vector{Float64} : optimized solution


"""
function testjl_rot(ig      :: Vector{Float64},
                    gamma   :: Vector{Float64},
                    ph_id   :: Int64,
                    gv, DB  )

    gv.maxeval   		= 1024;	
    gv.obj_tol   		= 1e-6;	

    # set the chemical potential Gamma

    # set the initial guess (xeos)
    SS_ref_db   = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

    unsafe_copyto!(gv.gam_tot, pointer(gamma), gv.len_ox)

    SS_ref_db[ph_id] = LibMAGEMin.rotate_hyperplane(gv,SS_ref_db[ph_id])
    SS_ref_db[ph_id] = LibMAGEMin.restrict_SS_HyperVolume(	gv,SS_ref_db[ph_id],1.0	);

    unsafe_copyto!(SS_ref_db[ph_id].iguess, pointer(ig), SS_ref_db[ph_id].n_xeos)

    if ph_id == 1
    SS_ref_db[ph_id] = LibMAGEMin.NLopt_opt_ig_spn_function(gv,SS_ref_db[ph_id])
    elseif ph_id == 4
    SS_ref_db[ph_id] = LibMAGEMin.NLopt_opt_ig_cpx_function(gv,SS_ref_db[ph_id])
    elseif ph_id == 7
    SS_ref_db[ph_id] = LibMAGEMin.NLopt_opt_ig_hb_function(gv,SS_ref_db[ph_id])
    end

    x   = unsafe_wrap(Vector{Float64},SS_ref_db[ph_id].xeos,SS_ref_db[ph_id].n_xeos)

    return xsol = copy(x)
end

"""
    select_ph(ph)

    Select the phase to minimize

    # Arguments
    - ph :: String : phase name

    # Returns
    - gamma :: Vector{Float64} : initial chemical potentials
    - pc :: Vector{Float64} : initial guess
    - P :: Float64 : pressure
    - T :: Float64 : temperature
    - ph_id :: Int64 : phase id
"""
function select_ph(ph :: String)

    # 1 - spn, 4 - cpx, 7 - hb

    gamma       = Vector{Float64}(undef, 11)
    P, T, ph_id = Float64, Float64, Int64

    if ph == "spn"
        include("/home/seph42/julia_scripts/XMin/spn_pc_list.jl")
        gamma       = [-1011.909631, -1829.092564, -819.264126, -695.467358, -412.948568, -971.890270, -876.544354, -1073.640927, -276.590707, -1380.299631, 0.0]
        P           =  12.0
        T           =  1200.0
        ph_id       =  1
    end
    if ph == "cpx"
        include("/home/seph42/julia_scripts/XMin/cpx_pc_list.jl")
        gamma       = [-1011.909631, -1829.092564, -819.264126, -695.467358, -412.948568, -971.890270, -876.544354, -1073.640927, -276.590707, -1380.299631, 0.0]
        P           =  12.0
        T           =  1100.0
        ph_id       =  4
    end
    if ph == "hb"
        include("/home/seph42/julia_scripts/XMin/hb_pc_list.jl")
        gamma       = [-960.9655,    -1768.2476,   -788.4474,    -678.9683,    -355.2975,    -914.9708,    -839.9561,    -1008.3630,   -263.7269,    -1262.6087,   -368.4674]
        P           =  5.0
        T           =  650.0
        ph_id       =  7
    end
    
    pc          = pc_list();
    return gamma, pc, P, T, ph_id
end


# Initialize database  - new way
using MAGEMin_C, BenchmarkTools

data        =   Initialize_MAGEMin("ig", verbose=true);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);

ph          = "spn"
gamma, pc, P, T, ph_id = select_ph(ph)
out         =   point_wise_minimization(P,T, data);

if ph      == "hb"
    ig        = [0.6204951407643728, 0.33881752819644395, 0.04068733103918367, 0.6305632282659267, 0.36943677173407297, 0.4382921274417739, 0.14257898878269964, 0.27831816917201424, 0.10980305743829911, 0.031007657165212758, 0.9184137906311484, 0.009163076636147906, 0.046488434108646656, 0.025934698624057062, 0.7085266926153586, 0.2914733073846413, 0.9689923428347873];
    gamma_eq  = [-960.9655,    -1768.2476,   -788.4474,    -678.9683,    -355.2975,    -914.9708,    -839.9561,    -1008.3630,   -263.7269,    -1262.6087,   -368.4674]
    sol       = [0.3478656659211239, 0.27831800067674756, 0.025934664164169143, 0.37950481793022106, 0.10721169904978942, 0.9184139366839471, 0.10980310198258116, 0.03100769194612725, -0.02157122090105998, 0.1024081525008751]

elseif ph  == "spn"
    ig        = [0.6214899217498702, 0.10378242277604803, 0.25186128463101126, 0.022866370843070603, 0.11796061376230191, 0.02179983119726647, 0.8108683444311726, 0.0014850439357559522, 0.045489549450975345, 0.0023966172225273236]
    gamma_eq  = [-1011.909631, -1829.092564, -819.264126, -695.467358, -412.948568, -971.890270, -876.544354, -1073.640927, -276.590707, -1380.299631, 0.0]
    sol       = [0.2642958560657611, 0.03129706158032523, 0.09019727879183574, 0.012965436292352656, 0.4490739361546607, 0.11764764891307426, 0.038932447396348484]
elseif ph  == "cpx"
    ig        = [0.744974756737465, 0.054254956805322656, 0.14567733693852883, 0.02154640055829482, 0.015785376254575986, 0.01776117270581313, 0.16399672343680996, 0.04134566405506293, 0.668655157506766, 0.12062805970438997, 0.00537439529697088, 0.9537354979191687, 0.04626450208083163]
    gamma_eq  = [-1011.909631, -1829.092564, -819.264126, -695.467358, -412.948568, -971.890270, -876.544354, -1073.640927, -276.590707, -1380.299631, 0.0]
    sol       = [0.0951654952819037, 0.09252901937086032, 0.20534220739762649, 0.12062808412939652, -0.027281442080958698, 0.021546414721478037, 0.015785385179589818, 0.017761159878362867, 0.005374395625155299]
end


np      = 10000                 # number of tested minimization
tms = []
xeostot = []
gam = []
# @benchmark begin
    for i = 1:np
        gamma .= gamma_eq .+ 40.0*(rand()-0.5)
        push!(gam,Array(gamma));
        t = @elapsed (xeos = testjl_rot(ig, gamma, ph_id, data.gv[1], data.DB[1]))
        push!(tms,t)
        push!(xeostot,xeos)
    end
# end

plot()
d = zeros(np)
for i=1:np
    # d[i] = norm(sf_out[:,i] .- sf_in[:,i])
    d[i] = norm(gam[i] .- gamma_eq)
end
leg = string(np)*" points";
scatter!(d,tms*1e6,ma=0.4,mc=:gray,ms=2,xlabel="Norm(ΔΓ)",ylabel="Minimization time [µs]",title="spinel",legend=false, label=leg, dpi=300)
plot!(legend=:bottomright, legendcolumns=1)
ylims!((0.0,500))
savefig("tms_vs_normdG_"*ph*"_MAGEMin.png")


###################################################################

tms = []
xeostot = []
np =size(pc,1)
nt = 0
# @benchmark begin
    x_in   = Vector{Float64}[]
    x_out  = Vector{Float64}[]
    for i = 1:np
        # print("$i ")
        t = @elapsed (xeos = testjl(pc[i,:], gamma, ph_id, data.gv[1], data.DB[1]))
        
        if abs(norm(xeos) - norm(sol)) < 1e-4
            push!(tms,t)
            push!(x_in,pc[i,:])
            push!(x_out,xeos)
            nt  += 1
        end
        # print("$xeos\n")
    end
# end

print("std(tms) =  $(std(tms))\nmean(tms) =  $(mean(tms))\nminimum(tms) =  $(minimum(tms))\nmaximum(tms) =  $(maximum(tms))\n")




plot()
d = zeros(nt);
for i=1:nt
    d[i] = norm(x_out[i,:] .- x_in[i,:])
end
leg = string(nt)*" points";
scatter!(d,tms*1e6,ma=0.4,mc=:gray,ms=2,xlabel="Norm(SFᶠ - SFⁱ)",ylabel="Minimization time [µs]",title=ph,legend=false, label=leg, dpi=300)
plot!(legend=:bottomright, legendcolumns=1)
if ph == "hb"
    ylims!((150.0,1000))
elseif ph == "spn"
    ylims!((50.0,300))
elseif ph == "cpx"
    ylims!((50.0,800))
end

savefig("tms_vs_distance_"*ph*".png")



# test
# xsol = testjl(pc[1,:], gamma, ph_id, data.gv[1], data.DB[1])


# SLSQP SOL C
# Time  (mean ± σ):   82.781 μs ±   5.835 μs  ┊ GC (mean ± σ):  0.00% ± 0.00% #TOL 1e-6
# Time  (mean ± σ):   93.569 μs ±   5.829 μs  ┊ GC (mean ± σ):  0.00% ± 0.00% #TOL 1e-13

# NULL SOL JULIA
# Time  (mean ± σ):   89.116 μs ±   3.359 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%


# gv = LibMAGEMin.LP_jl(z_b,gv,PC_read,SS_objective,NLopt_opt,SS_ref_db,ph_id-1)
# SS_ref_db[ph_id] = LibMAGEMin.PC_function(gv, PC_read, SS_ref_db[ph_id], z_b, ph_id-1)
# x   = unsafe_wrap(Vector{Float64},SS_ref_db[ph_id].xeos,SS_ref_db[ph_id].n_xeos);
