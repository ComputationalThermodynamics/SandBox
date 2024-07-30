using PlotlyJS

include("values_plot.jl")

function generate_min_time_vs_distance_figure(ph, gv, gm)
    # Initialize variables
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
            ph.ig .= get_ig(ph.pc_list[ig, :])
            push!(sf_in, Array(ph.ig))
            t = @elapsed null_min!(gv, ph, gm)

            push!(sf_out, Array(gm.x))
            append!(minG, ph.G)
            append!(tms, t)
            append!(n_out, gm.o)
            append!(n_in, gm.i)
            append!(n_ite, gm.ite)
        end
    end

    sf_out = reduce(hcat, sf_out)
    sf_in  = reduce(hcat, sf_in)

    Plots.plot()
    d = zeros(np)
    for i = 1:np
        ph.sf .= sf_out[:, i]
        get_cv!(ph)
        cv1 = copy(ph.cv)
        ph.sf .= sf_in[:, i]
        get_cv!(ph)
        cv2 = copy(ph.cv)
        d[i] = norm(cv1 .- cv2)
    end

    println("std(tms) = $(std(tms))")
    println("mean(tms) = $(mean(tms))")
    println("minimum(tms) = $(minimum(tms))")
    println("maximum(tms) = $(maximum(tms))")

    leg = string(np) * " points"
    scatter!(d, tms * 1e6, ma = 0.4, mc = :gray, ms = 2, xlabel = "Norm(SFᶠ - SFⁱ)", ylabel = "Minimization time [µs]", title = ph.ph, legend = false, label = leg, dpi = 300)
    plot!(legend = :bottomright, legendcolumns = 1)
    if ph.ph == "clino-amphibole"
        ylims!((150.0, 1000))
    elseif ph.ph == "spinel"
        ylims!((50.0, 300))
    elseif ph.ph == "clinopyroxene"
        ylims!((50.0, 800))
    end

    Plots.savefig("tms_vs_distance_" * ph.ph * ".png")

    Plots.plot()
    leg = string(np) * " points"
    scatter!(n_ite, tms * 1e6, ma = 0.4, mc = :gray, ms = 2, xlabel = "Number of iterations", ylabel = "Minimization time [µs]", title = ph.ph, legend = false, label = leg, dpi = 300)
    plot!(legend = :bottomright, legendcolumns = 1)
    if ph.ph == "clino-amphibole"
        ylims!((150.0, 600))
    elseif ph.ph == "spinel"
        ylims!((50.0, 300))
    elseif ph.ph == "clinopyroxene"
        ylims!((50.0, 400))
    end

    Plots.savefig("tms_vs_nite_" * ph.ph * ".png")
end


function generate_min_time_vs_normDeltaGamma_figure(ph, gv, gm)
    # Initialize variables
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
            ph.ig .= get_ig(ph.pc_list[ig, :])
            push!(sf_in, Array(ph.ig))
            t = @elapsed null_min!(gv, ph, gm)

            push!(sf_out, Array(gm.x))
            append!(minG, ph.G)
            append!(tms, t)
            append!(n_out, gm.o)
            append!(n_in, gm.i)
            append!(n_ite, gm.ite)
        end
    end

    sf_out = reduce(hcat, sf_out)
    sf_in  = reduce(hcat, sf_in)

    Plots.plot()
    d = zeros(np)
    for i = 1:np
        ph.sf .= sf_out[:, i]
        get_cv!(ph)
        cv1 = copy(ph.cv)
        ph.sf .= sf_in[:, i]
        get_cv!(ph)
        cv2 = copy(ph.cv)
        d[i] = norm(cv1 .- cv2)
    end

    println("std(tms) = $(std(tms))")
    println("mean(tms) = $(mean(tms))")
    println("minimum(tms) = $(minimum(tms))")
    println("maximum(tms) = $(maximum(tms))")

    leg = string(np) * " points"
    scatter!(d, tms * 1e6, ma = 0.4, mc = :gray, ms = 2, xlabel = "Norm(SFᶠ - SFⁱ)", ylabel = "Minimization time [µs]", title = ph.ph, legend = false, label = leg, dpi = 300)
    plot!(legend = :bottomright, legendcolumns = 1)
    if ph.ph == "clino-amphibole"
        ylims!((150.0, 1000))
    elseif ph.ph == "spinel"
        ylims!((50.0, 300))
    elseif ph.ph == "clinopyroxene"
        ylims!((50.0, 800))
    end

    Plots.savefig("tms_vs_distance_" * ph.ph * ".png")

    Plots.plot()
    leg = string(np) * " points"
    scatter!(n_ite, tms * 1e6, ma = 0.4, mc = :gray, ms = 2, xlabel = "Number of iterations", ylabel = "Minimization time [µs]", title = ph.ph, legend = false, label = leg, dpi = 300)
    plot!(legend = :bottomright, legendcolumns = 1)
    if ph.ph == "clino-amphibole"
        ylims!((150.0, 600))
    elseif ph.ph == "spinel"
        ylims!((50.0, 300))
    elseif ph.ph == "clinopyroxene"
        ylims!((50.0, 400))
    end

    Plots.savefig("tms_vs_nite_" * ph.ph * ".png")
end

function create_box_plot()


    function generate_box(t_hb, t_cpx, t_spn, name, color)
        x0 = Vector{String}(undef, length(t_hb))
        x0 .= "clino-amphibole"

        x1 = Vector{String}(undef, length(t_cpx))
        x1 .= "clinopyroxene"

        x2 = Vector{String}(undef, length(t_spn))
        x2 .= "spinel"

        x = [x0; x1; x2]

        box(y = [(t_hb .* 1e3); (t_cpx .* 1e3); (t_spn .* 1e3)],
            boxpoints = "all",
            x = x,
            name = name,
            marker_color = color,
            marker = PlotlyJS.attr(size = 1.0, color = color))
    end

    CCSAQ = generate_box(t_hb_ccsaq, t_cpx_ccsaq, t_spn_ccsaq, "CCSAQ", c[1])
    SLSQP = generate_box(t_hb_slsqp, t_cpx_slsqp, t_spn_slsqp, "SLSQP", c[2])
    CG = generate_box(t_hb_null_cg, t_cpx_null_cg, t_spn_null_cg, "CG", c[3])
    BFGS = generate_box(t_hb_null_bfgs, t_cpx_null_bfgs, t_spn_null_bfgs, "BFGS", c[4])

    layout = Layout(title = PlotlyJS.attr(
                        x = 0.5,
                        xanchor = "center",
                        yanchor = "top"
                    ),
                    yaxis_title = "Time [ms]",
                    yaxis_type = "log",
                    xaxis_title = "Solution model",
                    boxmode = "group",
                    yaxis_range = [-1.5, 1.0],
                    template = "plotly_white")

    figure = PlotlyJS.plot([CCSAQ, SLSQP, CG, BFGS], layout)

    PlotlyJS.savefig(figure, "minimization_time_plot.pdf")
end

function generate_figure_solvus(ph, gv, gm)
    if ph.ph == "spinel"
        minG = []
        tms = []
        n_out = []
        n_in = []
        n_ite = []
        sf_in = Vector{Float64}[]
        sf_out = Vector{Float64}[]
        trackSF = Vector{Float64}[]
        np = size(ph.pc_list)[1]

        Plots.plot()
        @time begin
            for ig in 1:np
                empty!(trackSF)
                ph.ig .= get_ig(ph.pc_list[ig, :])
                push!(sf_in, Array(ph.ig))
                t = @elapsed null_min_BFGS_s!(gv, ph, gm, trackSF)

                push!(sf_out, Array(gm.x))
                append!(minG, ph.G)
                append!(tms, t)
                append!(n_out, gm.o)
                append!(n_in, gm.i)
                append!(n_ite, gm.ite)

                line_x = zeros(length(trackSF))
                line_y = zeros(length(trackSF))
                for i in 1:length(trackSF)
                    line_x[i] = trackSF[i][1]
                    line_y[i] = trackSF[i][10]
                end

                if gm.x[10] > 0.2
                    plot!(line_x, line_y, lc=:deepskyblue3, lw=0.5, la=0.3, legend=false, dpi=600)
                else
                    plot!(line_x, line_y, lc=:gray, lw=0.5, la=0.3, legend=false, dpi=600)
                end
            end
        end

        sf_out = reduce(hcat, sf_out)
        sf_in = reduce(hcat, sf_in)

        scatter!(sf_in[1, :], sf_in[10, :], mc=:black, ms=1, xlabel="xMgT", ylabel="xTiM")
        scatter!(sf_out[1, :], sf_out[10, :], mc=:orange, markershape=:star, ms=5)

        Plots.savefig("solvus_" * ph.ph * ".png")
    end
end