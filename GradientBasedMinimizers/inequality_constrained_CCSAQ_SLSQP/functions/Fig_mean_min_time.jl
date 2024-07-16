
tms     = []
xeostot = []
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