# generate figure

np      = 10000                 # number of tested minimization
minG    = []
tms     = []
n_out   = []
n_in    = []
n_ite   = []
sf_in   = Vector{Float64}[]
sf_out  = Vector{Float64}[]
gam     = Vector{Float64}[]

# here we set the solution from minimization as starting point, and perturb it below
if join(ph.ph)      == "clino-amphibole"
    ph.ig        .= [0.6204951407643728, 0.33881752819644395, 0.04068733103918367, 0.6305632282659267, 0.36943677173407297, 0.4382921274417739, 0.14257898878269964, 0.27831816917201424, 0.10980305743829911, 0.031007657165212758, 0.9184137906311484, 0.009163076636147906, 0.046488434108646656, 0.025934698624057062, 0.7085266926153586, 0.2914733073846413, 0.9689923428347873];
    gamma_eq      = [-960.9655,    -1768.2476,   -788.4474,    -678.9683,    -355.2975,    -914.9708,    -839.9561,    -1008.3630,   -263.7269,    -1262.6087,   -368.4674]
elseif join(ph.ph)  == "spinel"
    ph.ig        .= [0.6214899217498702, 0.10378242277604803, 0.25186128463101126, 0.022866370843070603, 0.11796061376230191, 0.02179983119726647, 0.8108683444311726, 0.0014850439357559522, 0.045489549450975345, 0.0023966172225273236]
    gamma_eq      = [-1011.909631, -1829.092564, -819.264126, -695.467358, -412.948568, -971.890270, -876.544354, -1073.640927, -276.590707, -1380.299631, 0.0]
elseif join(ph.ph)  == "clinopyroxene"
    ph.ig        .= [0.744974756737465, 0.054254956805322656, 0.14567733693852883, 0.02154640055829482, 0.015785376254575986, 0.01776117270581313, 0.16399672343680996, 0.04134566405506293, 0.668655157506766, 0.12062805970438997, 0.00537439529697088, 0.9537354979191687, 0.04626450208083163]
    gamma_eq      = [-1011.909631, -1829.092564, -819.264126, -695.467358, -412.948568, -971.890270, -876.544354, -1073.640927, -276.590707, -1380.299631, 0.0]
end

@time begin
    for ig = 1:np
        ph.gamma .= gamma_eq .+ 40.0*(rand()-0.5)
        get_gb!(ph);

        push!(sf_in,Array(ph.ig));
        t = @elapsed null_min!(gv,ph,gm);
        push!(gam,Array(ph.gamma));
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
d = zeros(np)
for i=1:np
    # d[i] = norm(sf_out[:,i] .- sf_in[:,i])
    d[i] = norm(gam[i] .- gamma_eq)
end
leg = string(np)*" points";
scatter!(d,tms*1e6,ma=0.4,mc=:gray,ms=2,xlabel="Norm(ΔΓ)",ylabel="Minimization time [µs]",title=join(ph.ph),legend=false, label=leg, dpi=300)
plot!(legend=:bottomright, legendcolumns=1)

if join(ph.ph) == "clino-amphibole"
    ylims!((0.0,600))
elseif join(ph.ph) == "spinel"
    ylims!((25.0,200))
elseif join(ph.ph) == "clinopyroxene"
    ylims!((50.0,400))
end

savefig("tms_vs_normdG_"*join(ph.ph)*".png")
