
"""
    Fig_mean_min_time(ph, gamma, ph_id, ig, gamma_eq)

Generate a scatter plot of minimization time versus the norm of ΔΓ for a given phase.
"""
function Fig_mean_min_time(ph, gamma, ph_id, ig, gamma_eq)

    np      = 10000
    tms     = []
    xeostot = []
    gam     = []

    for _ = 1:np
        gamma .= gamma_eq .+ 40.0*(rand()-0.5)
        push!(gam,Array(gamma));
        t = @elapsed (xeos = testjl_rot(ig, gamma, ph_id, data.gv[1], data.DB[1]))
        push!(tms,t)
        push!(xeostot,xeos)
    end


    plot()
    d = zeros(np)
    for i=1:np
        d[i] = norm(gam[i] .- gamma_eq)
    end
    leg = string(np)*" points";
    scatter!(d,tms*1e6,ma=0.4,mc=:gray,ms=2,xlabel="Norm(ΔΓ)",ylabel="Minimization time [µs]",title="$(String(ph))",legend=false, label=leg, dpi=300)
    plot!(legend=:bottomright, legendcolumns=1)
    ylims!((0.0,500))


    savefig("tms_vs_normdG_"*String(ph)*"_MAGEMin.png")

    return nothing
end

"""
    Fig_minTime_vs_normDeltaGamma(ph, gamma, pc, ph_id, ig, sol)

Generate a scatter plot of minimization time versus the norm (SFᶠ - SFⁱ) for a given phase.
"""
function Fig_minTime_vs_normDeltaGamma(ph, gamma, pc, ph_id, ig, sol)

    tms     = []
    np      = size(pc,1)
    nt      = 0

    x_in   = Vector{Float64}[]
    x_out  = Vector{Float64}[]
    for i = 1:np
        t = @elapsed (xeos = testjl(pc[i,:], gamma, ph_id, data.gv[1], data.DB[1]))

        if abs(norm(xeos) - norm(sol)) < 1e-4
            push!(tms,t)
            push!(x_in,pc[i,:])
            push!(x_out,xeos)
            nt  += 1
        end
    end

    print("std(tms) =  $(std(tms))\nmean(tms) =  $(mean(tms))\nminimum(tms) =  $(minimum(tms))\nmaximum(tms) =  $(maximum(tms))\n")

    plot()
    d = zeros(nt);
    for i=1:nt
        d[i] = norm(x_out[i,:] .- x_in[i,:])
    end
    leg = string(nt)*" points";
    scatter!(d,tms*1e6,ma=0.4,mc=:gray,ms=2,xlabel="Norm(SFᶠ - SFⁱ)",ylabel="Minimization time [µs]",title=String(ph),legend=false, label=leg, dpi=300)
    plot!(legend=:bottomright, legendcolumns=1)
    if ph == :hb
        ylims!((150.0,1000))
    elseif ph == :spn
        ylims!((50.0,300))
    elseif ph == :cpx
        ylims!((50.0,800))
    end

    savefig("tms_vs_distance_"*String(ph)*".png")

    return nothing
end