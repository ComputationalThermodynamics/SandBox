# generate figure

minG    = []
tms     = []
n_out   = []
n_in    = []
n_ite   = []
sf_in   = Vector{Float64}[]
sf_out  = Vector{Float64}[]
np      = size(ph.pc_list)[1]

@time begin

    for ig = 1:np
        ph.ig          .= get_ig(ph.pc_list[ig,:]);
        push!(sf_in,Array(ph.ig));
        t = @elapsed null_min!(gv,ph,gm);

        push!(sf_out,Array(gm.x));
        append!(minG,ph.G);
        append!(tms,t);
        append!(n_out,gm.o);
        append!(n_in,gm.i);
        append!(n_ite,gm.ite);

    end
end
sf_out = reduce(hcat,sf_out);
sf_in  = reduce(hcat,sf_in);


plot()
d = zeros(np);
for i=1:np

    ph.sf .= sf_out[:,i]
    get_cv!(ph)
    cv1 = copy(ph.cv)
    ph.sf .= sf_in[:,i]
    get_cv!(ph)
    cv2 = copy(ph.cv)
    d[i] = norm(cv1 .- cv2)

    # d[i] = norm(sf_out[:,i] .- sf_in[:,i])

end

print("std(tms) =  $(std(tms))\nmean(tms) =  $(mean(tms))\nminimum(tms) =  $(minimum(tms))\nmaximum(tms) =  $(maximum(tms))\n")


leg = string(np)*" points";
scatter!(d,tms*1e6,ma=0.4,mc=:gray,ms=2,xlabel="Norm(SFᶠ - SFⁱ)",ylabel="Minimization time [µs]",title=ph.ph,legend=false, label=leg, dpi=300)
plot!(legend=:bottomright, legendcolumns=1)
if ph.ph == "clino-amphibole"
    ylims!((150.0,1000))
elseif ph.ph == "spinel"
    ylims!((50.0,300))
elseif ph.ph == "clinopyroxene"
    ylims!((50.0,800))
end

savefig("tms_vs_distance_"*ph.ph*".png")

plot()
leg = string(np)*" points";
scatter!(n_ite,tms*1e6,ma=0.4,mc=:gray,ms=2,xlabel="Number of iterations",ylabel="Minimization time [µs]",title=ph.ph,legend=false, label=leg, dpi=300)
plot!(legend=:bottomright, legendcolumns=1)
if ph.ph == "clino-amphibole"
    ylims!((150.0,600))
elseif ph.ph == "spinel"
    ylims!((50.0,300))
elseif ph.ph == "clinopyroxene"
    ylims!((50.0,400))
end

savefig("tms_vs_nite_"*ph.ph*".png")
