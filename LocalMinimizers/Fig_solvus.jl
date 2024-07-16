# generate figure

if ph.ph == "spinel"
    minG    = []
    tms     = []
    n_out   = []
    n_in    = []
    n_ite   = []
    sf_in   = Vector{Float64}[]
    sf_out  = Vector{Float64}[]
    trackSF         = Vector{Float64}[];
    np      = size(ph.pc_list)[1]

    plot();
    @time begin

        for ig = 1:np

            empty!(trackSF)
            ph.ig          .= get_ig(ph.pc_list[ig,:]);
            push!(sf_in,Array(ph.ig));
            t = @elapsed null_min_BFGS_s!(gv,ph,gm,trackSF);

            push!(sf_out,Array(gm.x));
            append!(minG,ph.G);
            append!(tms,t);
            append!(n_out,gm.o);
            append!(n_in,gm.i);
            append!(n_ite,gm.ite);

            line_x = zeros(length(trackSF))
            line_y = zeros(length(trackSF))
            for i=1:length(trackSF)
                line_x[i] = trackSF[i,][1]
                line_y[i] = trackSF[i,][10]
            end
            # line_x = [sf_in[1,i],sf_out[1,i]]
            # line_y = [sf_in[10,i],sf_out[10,i]]
            if gm.x[10] > 0.2
                plot!(line_x,line_y,lc=:deepskyblue3,lw=0.5,la=0.3,legend=false,dpi=600)
            else
                plot!(line_x,line_y,lc=:gray,lw=0.5,la=0.3,legend=false,dpi=600)
            end

        end
    end
    sf_out = reduce(hcat,sf_out);
    sf_in  = reduce(hcat,sf_in);

    scatter!(sf_in[1,:],sf_in[10,:],mc=:black,ms=1,xlabel="xMgT",ylabel="xTiM")
    scatter!(sf_out[1,:],sf_out[10,:],mc=:orange,markershape=:star,ms=5)

    savefig("solvus_"*ph.ph*".png")
end
