
np      = 10000 
tms     = []
xeostot = []
gam     = []

for i = 1:np
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
scatter!(d,tms*1e6,ma=0.4,mc=:gray,ms=2,xlabel="Norm(ΔΓ)",ylabel="Minimization time [µs]",title="spinel",legend=false, label=leg, dpi=300)
plot!(legend=:bottomright, legendcolumns=1)
ylims!((0.0,500))


savefig("tms_vs_normdG_"*ph*"_MAGEMin.png")
