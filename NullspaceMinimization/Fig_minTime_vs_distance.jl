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
    d[i] = norm(sf_out[:,i] .- sf_in[:,i])
end
leg = string(np)*" points";
scatter!(d,tms*1e6,ma=0.4,mc=:gray,ms=2,xlabel="Norm(SFᶠ - SFⁱ)",ylabel="Minimization time [µs]",title=join(ph.ph),legend=false, label=leg, dpi=300)
plot!(legend=:bottomright, legendcolumns=1)
if join(ph.ph) == "clino-amphibole"
    ylims!((150.0,500))
elseif join(ph.ph) == "spinel"
    ylims!((50.0,300))
elseif join(ph.ph) == "clinopyroxene"
    ylims!((100.0,500))
end

savefig("tms_vs_distance_"*join(ph.ph)*".png")
