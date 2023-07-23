# generate figure

phase       = ["clino-amphibole","clinopyroxene","spinel"]
mean_minG   = [-1.3191066463022626e-5,3.5850555168642084e-7,-1.742843796266652e-5]
std_minG    = [7.561616392069336e-14,1.1880974573794076e-15,1.0355434263191188e-13]

mean_tms    = [0.0003167131692929292,0.0002539734131910765,0.00014003279618671934]      #sec
std_tms     = [0.0001404211479561104,0.0001471354718169665,3.9017098592550434e-5]       #sec


plot()
counts = mean_tms.*1e6;
barswitherror = bar(phase, counts, legend=false,label=false,alpha=0.0);
leg = "Nullspace";
scatter!(phase, counts, yerr=std_tms.*1e6, ylabel="minimization time [µs]", title="CCSAQ vs Nullspace", marker=stroke(1.5), legend=false, label=leg, dpi=300)
plot!(legend=:bottomright, legendcolumns=1)
ylims!((0.0,500))
savefig("CCSAQ_vs_NLopt_min_time.png")
