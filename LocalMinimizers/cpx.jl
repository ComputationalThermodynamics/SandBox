# Clinopyroxene data (Igneous database, Holland et al., 2018)
using StaticArrays;
using LinearAlgebra;
using BenchmarkTools


include("cpx_pc_list.jl")

"""
    Structure to store solution phase data
"""
struct solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W,  
                    n_ox_t_n_em, n_em_t_n_em, n_sf_t_n_sf, n_eq_t_n_sf, n_em_t_n_W, n_em_t_n_sf,    # number of entries in matrixes
                    N, _C, _T, _I}
    ph          ::  String       # name of phase
    n_eq_off    ::  MVector{n_sf, _I}
    P           ::  _T
    T           ::  _T

    bulk        ::  MVector{n_ox, _T} 
    gb          ::  MVector{n_em, _T}
    g0          ::  MVector{n_em, _T}
    v_nem       ::  MVector{n_em, _T}
    v_nsf       ::  MVector{n_sf, _T}
    v_E         ::  MMatrix{n_em, n_em, _T, n_em_t_n_em}
    v_A         ::  MMatrix{n_sf, n_sf, _T, n_sf_t_n_sf}
    gamma       ::  MVector{n_ox, _T}  

    emC         ::  MMatrix{n_em, n_ox, _T, n_ox_t_n_em}
    A           ::  MMatrix{n_eq, n_sf, _T, n_eq_t_n_sf}
    v           ::  MVector{n_em,_T}
    W           ::  MVector{n_W,_T}
    v_W         ::  MMatrix{n_em,n_W,_T, n_em_t_n_W} 

    ig          ::  MVector{n_sf,_T}
    sf          ::  MVector{n_sf,_T}
    cv          ::  MVector{n_xeos,_T}
    p           ::  MVector{n_em,_T}
    f           ::  MVector{1,_T}
    df          ::  MVector{n_em,_T}
    mat_phi     ::  MVector{n_em,_T}
    Gex         ::  MVector{n_em,_T}
    idm         ::  MVector{n_em,_T}
    
    mu          ::  MVector{n_em,_T}
    Graw        ::  MVector{1,_T}
    G           ::  MVector{1,_T}

    dpdsf       ::  MMatrix{n_em,n_sf, _T, n_em_t_n_sf}
    dGdsf       ::  MVector{n_sf, _T}
    grad        ::  MVector{n_sf, _T}
    sf_off      ::  MVector{n_sf, _T}
    pc_list     ::  Matrix{_T}

end


"""
    Get constant parameters
"""
function init_phase(gam,test)

    ph          = "clinopyroxene"; 
    n_ox        = 11
    n_sf        = 13;
    n_eq        = 4;
    n_eq_off    = MVector{n_sf}(zeros(Int64,n_sf));
    n_em        = 10;
    n_xeos      = 9;
    n_W         = 45;   # perhaps this can be computed somehow from the variables above?

    P           = 12.0;                                  # pressure kbar
    T           = 1100.0 + 273.15;                       # temperature C

    bulk        = MVector{n_ox}([0.38494,    0.01776,    0.02824,    0.50566,    0.05886,    0.0001,   0.00250,    0.0010,    0.00096,    0.00109,    0.0]);
    bulk       .= bulk./(sum(bulk));

    # reference Gibbs energy of em, and Gamma at which the equilibrium was computed (obtained using MAGEMin)
    g0          = SVector{n_em}([-3532.74915,  -2793.12846,  -3635.49886,  -3384.95041,  -3250.67812,  -3606.43710,  -3345.42582,  -3408.36774,  -3105.14810,  -3360.74459]);
    gamma_eq    = SVector{n_ox}([-1011.909631, -1829.092564, -819.264126, -695.467358, -412.948568, -971.890270, -876.544354, -1073.640927, -276.590707, -1380.299631, 0.0]);
    gamma_em    = SVector{n_ox}([-1016.399486, -1823.923434, -805.068209, -692.021970, -383.138855, -897.240517, -828.503925, -1062.477090, -310.231496, -1341.449679, 0.0]);
    gamma_pc    = SVector{n_ox}([-1006.885527, -1833.130473, -821.014905, -697.963191, -415.466554, -947.276954, -901.557511, -1058.244567, -236.176855, -1359.230803, 0.0]);

    emC         = SMatrix{n_ox, n_em}(
                   [    +2.000000  +0.000000  +1.000000  +1.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000;
                        +2.000000  +0.000000  +0.000000  +0.000000  +2.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000;
                        +1.000000  +1.000000  +1.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000;
                        +1.000000  +0.500000  +1.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.500000  +0.000000;
                        +1.000000  +0.500000  +1.000000  +0.000000  +1.000000  +0.000000  +0.000000  +0.000000  +0.500000  +0.000000  +0.000000;
                        +1.000000  +0.500000  +1.000000  +0.500000  +0.000000  +0.000000  +0.000000  +0.500000  +0.000000  +0.000000  +0.000000;
                        +2.000000  +0.500000  +0.000000  +0.000000  +0.000000  +0.000000  +0.500000  +0.000000  +0.000000  +0.000000  +0.000000;
                        +2.000000  +0.000000  +0.000000  +2.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000;
                        +2.000000  +0.000000  +0.000000  +1.000000  +1.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000;
                        +2.000000  +0.500000  +0.000000  +0.000000  +0.000000  +0.500000  +0.000000  +0.000000  +0.000000  +0.000000  +0.000000 ]
                        );


    A           = SMatrix{n_eq, n_sf}(
                    [ 1.0  1    1    1    1 1    0    0    0 0 0 0 0;
                    0  0    0    0    0 0    1    1    1 1 1 0 0;
                    -1 -1 -1/2 -1/2 -1/2 0  1/2  1/2  1/2 0 0 1 0;
                    1  1  1/2  1/2  1/2 0 -1/2 -1/2 -1/2 0 0 0 1]
                    );

    v            = SVector{n_em}([1.2, 1.0, 1.9, 1.9, 1.9, 1.9, 1.2, 1.0, 1.0, 1.2]);
    W            = SVector{n_W}([25.8, 13.0 - 0.06*P, 8.0, 8.0, 8.0, 26.0, 29.8, 20.6, 26.0, 25.0 - 0.1*P, 38.3, 43.3, 24.0, 24.0, 2.3, 3.5, 24.0, 2.0, 2.0, 6.0, 6.0, 45.2 - 0.35*P, 27.0 - 0.1*P, 6.0, 2.0, 6.0, 3.0, 52.3, 40.3, 3.0, 6.0, 3.0, 57.3, 45.3, 3.0, 16.0, 24.0, 22.0, 16.0, 40.0, 40.0, 28.0, 4.0, 40.0, 40.0]);
    v_W          = MMatrix{n_em,n_W}(zeros(n_em,n_W));
    
    # get partial derivatives of endmember fraction as function of site fraction (constant for amphibole)
    dpdsf       = MMatrix{n_em,n_sf}(zeros(n_em,n_sf));

    df          = MVector{n_em}(zeros(n_em));

    f           = MVector{1}(zeros(1))
    Graw        = MVector{1}(zeros(1))
    G           = MVector{1}(zeros(1))
    
    v_nem       = MVector{n_em}(zeros(n_em));
    v_nsf       = MVector{n_sf}(zeros(n_sf));
    v_E         = SMatrix{n_em, n_em}(Matrix(1.0I, n_em, n_em));
    v_A         = MMatrix{n_sf, n_sf}(zeros(n_sf,n_sf));

    gamma       = MVector{n_ox}(zeros(n_ox));
    gb          = MVector{n_em}(zeros(n_em));
    mat_phi     = MVector{n_em}(zeros(n_em));
    Gex         = MVector{n_em}(zeros(n_em));
    idm         = MVector{n_em}(zeros(n_em));
    mu          = MVector{n_em}(zeros(n_em));
    p           = MVector{n_em}(zeros(n_em));

    ig          = MVector{n_sf}(zeros(n_sf));
    sf          = MVector{n_sf}(zeros(n_sf));
    cv          = MVector{n_xeos}(zeros(n_xeos));
    dGdsf       = MVector{n_sf}(zeros(n_sf));
    grad        = MVector{n_sf}(zeros(n_sf));
    sf_off      = MVector{n_sf}(zeros(n_sf));

    if (gam == "em")
        gamma .= gamma_em;
    elseif  (gam == "pc")
        gamma .= gamma_pc;
    elseif  (gam == "eq")
        gamma .= gamma_eq;
    end

    pcs         = pc_list()         # appears slow (type instable?)

    cpx =  solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W, 
                         n_ox*n_em, n_em*n_em, n_sf*n_sf, n_eq*n_sf, n_em*n_W, n_em*n_sf,  
                         length(ph), Char, Float64, Int64}(   
                            ph,
                            n_eq_off,
                            
                            P,
                            T,
                            bulk, 
                            gb,
                            g0,

                            v_nem,
                            v_nsf,
                            v_E,
                            v_A,

                            gamma,
                            emC,
                            A,
                            v,
                            W,
                            v_W,

                            ig,
                            sf,
                            cv,
                            p,
                            f,
                            df,
                            mat_phi,
                            Gex,
                            idm,
                            mu,
                            Graw,
                            G,
                            dpdsf,
                            dGdsf,
                            grad,
                            sf_off,
                            pcs    );

    # change of base with respect to the Gibbs hyperplane given
    # this function update g0 and the subsequent minimization is achieved with respect to the Gibbs hyperplane given by gamma
    get_gb!(cpx);
    get_v_W!(cpx);

    return cpx;
end

"""
    calculate normalizing factor and needed partial derivatives
"""
function get_v_W!(ph::solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W}) where {n_ox, n_sf, n_eq, n_em, n_xeos, n_W}

    v_W = zeros(Float64, n_em, n_W)
    for i=1:n_em #lastindex(ph.p)
        it = 1;
        for j=1:n_em-1
            for k=j+1:n_em
                v_W[i,it] = (ph.W[it]*2.0*ph.v[i]/(ph.v[j]+ph.v[k]));
                it = it + 1;
            end
        end
    end
    ph.v_W .= v_W


    return nothing
end

"""
    calculate normalizing factor and needed partial derivatives
"""
function get_f!(ph,gv)

    #mul!(ph.v_nem,ph.emC,gv.apo);
    ph.v_nem .= ph.emC*gv.apo

    v         = ph.p'*ph.v_nem
    ph.f     .= (ph.bulk'*gv.apo)/(v);
    ph.df    .= (ph.v_nem)./(v);

    return nothing
end


"""
    Computes endmember proportions
"""
function get_p!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}

    for i=1:lastindex(ph.sf)
        if (ph.sf_off == 1)
            ph.sf[i] = 0.0;
        end
    end
    
    ph.p[1]          = -2*ph.sf[13] - ph.sf[8] - ph.sf[11] - ph.sf[7] - ph.sf[10] + 1;
    ph.p[2]          = -2*ph.sf[13]*(ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - 2*ph.sf[13]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) + ph.sf[2]/(ph.sf[2] + ph.sf[1]) - ph.sf[11]*(ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - ph.sf[11]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) - ph.sf[10]*(ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - ph.sf[10]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) + ph.sf[6]*(ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + ph.sf[6]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + (ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]);
    ph.p[3]          = 2*ph.sf[13] - ph.sf[5] - ph.sf[4] - 2*ph.sf[6];
    ph.p[4]          = ph.sf[5];
    ph.p[5]          = ph.sf[4];
    ph.p[6]          = 2*ph.sf[6];
    ph.p[7]          = ph.sf[10];
    ph.p[8]          = -2*ph.sf[13]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) + ph.sf[2]/(ph.sf[2] + ph.sf[1]) + ph.sf[8] - ph.sf[11]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) + ph.sf[7] - ph.sf[10]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) + ph.sf[6]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]);
    ph.p[9]          = 2*ph.sf[13]*(ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + 4*ph.sf[13]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) - 2*ph.sf[2]/(ph.sf[2] + ph.sf[1]) + ph.sf[11]*(ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + 2*ph.sf[11]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) + ph.sf[10]*(ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + 2*ph.sf[10]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) - ph.sf[6]*(ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - 2*ph.sf[6]*(ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])) - 2*(-ph.sf[2] - ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - (ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]);
    ph.p[10]          = ph.sf[11];


    return nothing
end


"""
    Computes partial derivative of endmember fraction as function of site fraction
"""
function get_dpdsf!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}

    ph.dpdsf[1,1] = 0;      ph.dpdsf[1,2] = 0;      ph.dpdsf[1,3] = 0;      ph.dpdsf[1,4] = 0;      ph.dpdsf[1,5] = 0;      ph.dpdsf[1,6] = 0;      ph.dpdsf[1,7] = -1;      ph.dpdsf[1,8] = -1;      ph.dpdsf[1,9] = 0;      ph.dpdsf[1,10] = -1;      ph.dpdsf[1,11] = -1;      ph.dpdsf[1,12] = 0;      ph.dpdsf[1,13] = -2;      
    ph.dpdsf[2,1] = ph.sf[2]*(2*ph.sf[13] + ph.sf[11] + ph.sf[10] - ph.sf[6] - 1)/(ph.sf[2]^2 + 2*ph.sf[2]*ph.sf[1] + ph.sf[1]^2);      ph.dpdsf[2,2] = ph.sf[1]*(-2*ph.sf[13] - ph.sf[11] - ph.sf[10] + ph.sf[6] + 1)/(ph.sf[2]^2 + 2*ph.sf[2]*ph.sf[1] + ph.sf[1]^2);      ph.dpdsf[2,3] = 0;      ph.dpdsf[2,4] = 0;      ph.dpdsf[2,5] = 0;      ph.dpdsf[2,6] = ph.sf[2]/(ph.sf[2] + ph.sf[1]);      ph.dpdsf[2,7] = 0;      ph.dpdsf[2,8] = 0;      ph.dpdsf[2,9] = 0;      ph.dpdsf[2,10] = -ph.sf[2]/(ph.sf[2] + ph.sf[1]);      ph.dpdsf[2,11] = -ph.sf[2]/(ph.sf[2] + ph.sf[1]);      ph.dpdsf[2,12] = 0;      ph.dpdsf[2,13] = -2*ph.sf[2]/(ph.sf[2] + ph.sf[1]);      
    ph.dpdsf[3,1] = 0;      ph.dpdsf[3,2] = 0;      ph.dpdsf[3,3] = 0;      ph.dpdsf[3,4] = -1;      ph.dpdsf[3,5] = -1;      ph.dpdsf[3,6] = -2;      ph.dpdsf[3,7] = 0;      ph.dpdsf[3,8] = 0;      ph.dpdsf[3,9] = 0;      ph.dpdsf[3,10] = 0;      ph.dpdsf[3,11] = 0;      ph.dpdsf[3,12] = 0;      ph.dpdsf[3,13] = 2;      
    ph.dpdsf[4,1] = 0;      ph.dpdsf[4,2] = 0;      ph.dpdsf[4,3] = 0;      ph.dpdsf[4,4] = 0;      ph.dpdsf[4,5] = 1;      ph.dpdsf[4,6] = 0;      ph.dpdsf[4,7] = 0;      ph.dpdsf[4,8] = 0;      ph.dpdsf[4,9] = 0;      ph.dpdsf[4,10] = 0;      ph.dpdsf[4,11] = 0;      ph.dpdsf[4,12] = 0;      ph.dpdsf[4,13] = 0;      
    ph.dpdsf[5,1] = 0;      ph.dpdsf[5,2] = 0;      ph.dpdsf[5,3] = 0;      ph.dpdsf[5,4] = 1;      ph.dpdsf[5,5] = 0;      ph.dpdsf[5,6] = 0;      ph.dpdsf[5,7] = 0;      ph.dpdsf[5,8] = 0;      ph.dpdsf[5,9] = 0;      ph.dpdsf[5,10] = 0;      ph.dpdsf[5,11] = 0;      ph.dpdsf[5,12] = 0;      ph.dpdsf[5,13] = 0;      
    ph.dpdsf[6,1] = 0;      ph.dpdsf[6,2] = 0;      ph.dpdsf[6,3] = 0;      ph.dpdsf[6,4] = 0;      ph.dpdsf[6,5] = 0;      ph.dpdsf[6,6] = 2;      ph.dpdsf[6,7] = 0;      ph.dpdsf[6,8] = 0;      ph.dpdsf[6,9] = 0;      ph.dpdsf[6,10] = 0;      ph.dpdsf[6,11] = 0;      ph.dpdsf[6,12] = 0;      ph.dpdsf[6,13] = 0;      
    ph.dpdsf[7,1] = 0;      ph.dpdsf[7,2] = 0;      ph.dpdsf[7,3] = 0;      ph.dpdsf[7,4] = 0;      ph.dpdsf[7,5] = 0;      ph.dpdsf[7,6] = 0;      ph.dpdsf[7,7] = 0;      ph.dpdsf[7,8] = 0;      ph.dpdsf[7,9] = 0;      ph.dpdsf[7,10] = 1;      ph.dpdsf[7,11] = 0;      ph.dpdsf[7,12] = 0;      ph.dpdsf[7,13] = 0;      
    ph.dpdsf[8,1] = (-ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 + (ph.sf[2] + ph.sf[1])^2*(ph.sf[2] + ph.sf[8] + (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7])) + (ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 - (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[1])^2)*(2*ph.sf[13] + ph.sf[11] + ph.sf[10] - ph.sf[6]))/((ph.sf[2] + ph.sf[1])^2*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2);      ph.dpdsf[8,2] = (-ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 + (ph.sf[2] + ph.sf[1])^2*(ph.sf[2] + ph.sf[8] + (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7])) - (ph.sf[2] + ph.sf[1])^2*(ph.sf[8] + ph.sf[7] + 1)*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + (ph.sf[2] + ph.sf[1])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 + (2*ph.sf[13] + ph.sf[11] + ph.sf[10] - ph.sf[6])*(ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 - (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[1])^2 + (ph.sf[2] + ph.sf[1])^2*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - (ph.sf[2] + ph.sf[1])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2))/((ph.sf[2] + ph.sf[1])^2*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2);      ph.dpdsf[8,3] = 0;      ph.dpdsf[8,4] = 0;      ph.dpdsf[8,5] = 0;      ph.dpdsf[8,6] = ph.sf[2]/(ph.sf[2] + ph.sf[1]) - (ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]);      ph.dpdsf[8,7] = (-2*ph.sf[13]*(ph.sf[2] + ph.sf[8]) + ph.sf[2] + ph.sf[8] - ph.sf[11]*(ph.sf[2] + ph.sf[8]) - ph.sf[10]*(ph.sf[2] + ph.sf[8]) + ph.sf[6]*(ph.sf[2] + ph.sf[8]) + (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7]) - (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + (ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2)/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2;      ph.dpdsf[8,8] = (ph.sf[2] + ph.sf[8] + (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7]) - (ph.sf[1] + ph.sf[7])*(-2*ph.sf[13] - ph.sf[11] - ph.sf[10] + ph.sf[6]) + (ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 - (ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])*(ph.sf[2] + 2*ph.sf[8] + ph.sf[7] + 1))/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2;      ph.dpdsf[8,9] = 0;      ph.dpdsf[8,10] = -ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]);      ph.dpdsf[8,11] = -ph.sf[2]/(ph.sf[2] + ph.sf[1]) + (ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]);      ph.dpdsf[8,12] = 0;      ph.dpdsf[8,13] = 2*(-ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[1]))/((ph.sf[2] + ph.sf[1])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]));      
    ph.dpdsf[9,1] = (2*ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 + (ph.sf[2] + ph.sf[1])^2*(-2*ph.sf[13]*(ph.sf[2] + ph.sf[8]) - ph.sf[2] - ph.sf[8] - ph.sf[11]*(ph.sf[2] + ph.sf[8]) - ph.sf[10]*(ph.sf[2] + ph.sf[8]) + ph.sf[6]*(ph.sf[2] + ph.sf[8]) - (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7])) + 2*(ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 - (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[1])^2)*(-2*ph.sf[13] - ph.sf[11] - ph.sf[10] + ph.sf[6]))/((ph.sf[2] + ph.sf[1])^2*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2);      ph.dpdsf[9,2] = (2*ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 + (ph.sf[2] + ph.sf[1])^2*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])*(2*ph.sf[13] + ph.sf[8] + ph.sf[11] + ph.sf[7] + ph.sf[10] - ph.sf[6] + 1) + (ph.sf[2] + ph.sf[1])^2*(-2*ph.sf[13]*(ph.sf[2] + ph.sf[8]) - ph.sf[2] - ph.sf[8] - ph.sf[11]*(ph.sf[2] + ph.sf[8]) - ph.sf[10]*(ph.sf[2] + ph.sf[8]) + ph.sf[6]*(ph.sf[2] + ph.sf[8]) - (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7])) - 2*(ph.sf[2] + ph.sf[1])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 + 2*(-2*ph.sf[13] - ph.sf[11] - ph.sf[10] + ph.sf[6])*(ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2 - (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[1])^2 + (ph.sf[2] + ph.sf[1])^2*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - (ph.sf[2] + ph.sf[1])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2))/((ph.sf[2] + ph.sf[1])^2*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2);      ph.dpdsf[9,3] = 0;      ph.dpdsf[9,4] = 0;      ph.dpdsf[9,5] = 0;      ph.dpdsf[9,6] = (-2*ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) + (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[1]))/((ph.sf[2] + ph.sf[1])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]));      ph.dpdsf[9,7] = (2*ph.sf[13]*(ph.sf[2] + ph.sf[8]) - ph.sf[2] - ph.sf[8] + ph.sf[11]*(ph.sf[2] + ph.sf[8]) + ph.sf[10]*(ph.sf[2] + ph.sf[8]) - ph.sf[6]*(ph.sf[2] + ph.sf[8]) - (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7]) + (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]))/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2;      ph.dpdsf[9,8] = (-2*ph.sf[13]*(ph.sf[2] + ph.sf[8]) - ph.sf[2] - ph.sf[8] - ph.sf[11]*(ph.sf[2] + ph.sf[8]) - ph.sf[10]*(ph.sf[2] + ph.sf[8]) + ph.sf[6]*(ph.sf[2] + ph.sf[8]) - (ph.sf[2] + ph.sf[8])*(ph.sf[8] + ph.sf[7]) - 2*(ph.sf[1] + ph.sf[7])*(2*ph.sf[13] + ph.sf[11] + ph.sf[10] - ph.sf[6]) + (ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])*(2*ph.sf[13] + ph.sf[2] + 2*ph.sf[8] + ph.sf[11] + ph.sf[7] + ph.sf[10] - ph.sf[6] + 1))/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])^2;      ph.dpdsf[9,9] = 0;      ph.dpdsf[9,10] = (2*ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[1]))/((ph.sf[2] + ph.sf[1])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]));      ph.dpdsf[9,11] = (2*ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[1]))/((ph.sf[2] + ph.sf[1])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]));      ph.dpdsf[9,12] = 0;      ph.dpdsf[9,13] = 2*(2*ph.sf[2]*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]) - (ph.sf[2] + ph.sf[8])*(ph.sf[2] + ph.sf[1]))/((ph.sf[2] + ph.sf[1])*(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7]));      
    ph.dpdsf[10,1] = 0;      ph.dpdsf[10,2] = 0;      ph.dpdsf[10,3] = 0;      ph.dpdsf[10,4] = 0;      ph.dpdsf[10,5] = 0;      ph.dpdsf[10,6] = 0;      ph.dpdsf[10,7] = 0;      ph.dpdsf[10,8] = 0;      ph.dpdsf[10,9] = 0;      ph.dpdsf[10,10] = 0;      ph.dpdsf[10,11] = 1;      ph.dpdsf[10,12] = 0;      ph.dpdsf[10,13] = 0;      

    return nothing
end

"""
    Calculate G excess
"""
function get_Gex!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}

    ph.mat_phi  .= ph.p.*ph.v./(ph.p'*ph.v);

    a = ph.v_E;
    b = ph.mat_phi;
    c = ph.v_W;

    @inbounds for i=1:n_em
        it      = 1;
        Gex_v = 0.0;
        @inbounds for j=1:n_em-1
            tmp = (a[i,j]-b[j]);
            for k=j+1:n_em
                Gex_v -= tmp*(a[i,k]-b[k])*c[i,it];
                it = it + 1;
            end
        end
        ph.Gex[i]  = Gex_v;
    end

    return nothing
end

"""
    Computes ideal mixing activity
"""
function get_idm!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}

    for i=1:n_sf
        if (ph.sf_off[i] == 1.0)
            ph.sf[i] = 1.0;
        end
    end

    sf1 = sqrt(sqrt(ph.sf[12]));
    sf2 = sqrt(ph.sf[12]);
    sf3 = sqrt(sqrt(ph.sf[13])); 

    ph.idm[1]          = ph.sf[9]*ph.sf[1]*sf2;
    ph.idm[2]          = ph.sf[2]*ph.sf[8]*sf2;
    ph.idm[3]          = 1.4142*ph.sf[3]*sf3*ph.sf[9]*sf1;
    ph.idm[4]          = 1.4142*sf3*ph.sf[9]*ph.sf[5]*sf1;
    ph.idm[5]          = 1.4142*sf3*ph.sf[9]*ph.sf[4]*sf1;
    ph.idm[6]          = 2.8284*sf3*ph.sf[9]*sqrt(ph.sf[1])*sf1*sqrt(ph.sf[6]);
    ph.idm[7]          = ph.sf[3]*ph.sf[10]*sf2;
    ph.idm[8]          = ph.sf[1]*ph.sf[7]*sf2;
    ph.idm[9]          = ph.sf[8]*ph.sf[1]*sf2;
    ph.idm[10]         = ph.sf[3]*ph.sf[11]*sf2;

    return nothing
end


"""
    compute chemical potential of endmembers
"""
function get_mu!(ph::solution_phase{n_ox, n_sf, n_eq, n_em},gv) where {n_ox, n_sf, n_eq, n_em}
    for i=1:lastindex(ph.idm)
        if (ph.idm[i] == 0.0)
            ph.idm[i] = 1.0;
        end
    end

    for i=1:lastindex(ph.idm)
        ph.mu[i] =  gv.R*ph.T*real(log(ph.idm[i])) + ph.gb[i] + ph.Gex[i];
    end
    return nothing
end

"""
    Computes change of base for the Gibbs-hyperplane
"""
function get_gb!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}

    ph.v_nem    .= ph.emC*ph.gamma
    ph.gb       .= ph.g0 .- ph.v_nem;
    return nothing
end


function get_cv!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}
    # sf   = [xMgM1, xFeM1, xAlM1, xFe3M1, xCrM1, xTiM1, xMgM2, xFeM2, xCaM2, xNaM2, xKM2, xSiT, xAlT];
    ph.cv[1] = (ph.sf[2] + ph.sf[8])/(ph.sf[2] + ph.sf[8] + ph.sf[1] + ph.sf[7])
    ph.cv[2] = 2*ph.sf[13]
    ph.cv[3] = ph.sf[8] + ph.sf[7]
    ph.cv[4] = ph.sf[10]
    ph.cv[5] = -ph.cv[1] + ph.sf[2]/(ph.sf[2] + ph.sf[1])
    ph.cv[6] = ph.sf[4]
    ph.cv[7] = ph.sf[5]
    ph.cv[8] = ph.sf[6]
    ph.cv[9] = ph.sf[11]

    return nothing;
end



"""
    convert compositional variables to site fractions
    the original compositional variables names and site-fraction are used in order to make the code more readable
"""
function get_ig(cv)

    x = cv[1];
    y = cv[2];
    o = cv[3];
    n = cv[4];
    Q = cv[5];
    f = cv[6];
    cr = cv[7];
    t = cv[8];
    k = cv[9];
    
    xMgM1 = 1 - k - n - Q + t - x - y + k*Q + n*Q + (-Q)*t + k*x + n*x + (-t)*x + Q*y + x*y;
    xFeM1 = Q + x + (-k)*Q + (-n)*Q + Q*t + (-k)*x + (-n)*x + t*x + (-Q)*y + (-x)*y;
    xAlM1 = -cr - f + k + n + y - 2*t;
    xFe3M1 = f;
    xCrM1 = cr;
    xTiM1 = t;
    xMgM2 = o + Q + (-k)*Q + (-n)*Q + Q*t + (-o)*x + (-Q)*y;
    xFeM2 = -Q + k*Q + n*Q + (-Q)*t + o*x + Q*y;
    xCaM2 = 1 - k - n - o;
    xNaM2 = n;
    xKM2 = k;
    xSiT = 1 - 1/2*y;
    xAlT = 1/2*y;

    #site fractions
    sf   = [xMgM1, xFeM1, xAlM1, xFe3M1, xCrM1, xTiM1, xMgM2, xFeM2, xCaM2, xNaM2, xKM2, xSiT, xAlT];

    return sf
end


"""
    Computes Gibbs energy and first derivative of the solution phase
"""
function compute_G!(ph,gv)

    get_p!(ph);
    get_dpdsf!(ph);
    get_f!(ph,gv);
    get_Gex!(ph);
    get_idm!(ph);
    get_mu!(ph,gv);

    # Compute raw and normalized Gibbs energy
    ph.Graw   .= ph.mu'*ph.p;
    ph.G      .= ph.Graw[1]*ph.f[1];

    return nothing

end

function compute_G_dG!(gm,ph::solution_phase{n_ox, n_sf, n_eq, n_em},gv) where {n_ox, n_sf, n_eq, n_em}
    # total time of routine: 1.484 μs
    get_p!(ph);         # 23.012  ns
    get_dpdsf!(ph);     # 66.622  ns
    get_f!(ph,gv);      # 53.622  ns
    get_Gex!(ph);       # 514.728 ns -> 408.673 ns
    get_idm!(ph);       # 448.106 ns -> 136.606 ns -> 18.877 ns (using sqrt)
    get_mu!(ph,gv);     # 219.544 ns

    #   rest of routine: ~300 ns

    # Compute raw and normalized Gibbs energy
    ph.Graw   .= ph.mu'*ph.p;
    ph.G      .= ph.Graw[1]*ph.f[1];

    # compute first derivatives
    ph.v_nem  .= (ph.mu .- ph.df .*ph.Graw).*ph.f;
    ph.dGdsf  .= ph.dpdsf'*ph.v_nem
    #mul!(ph.dGdsf,ph.dpdsf',ph.v_nem);
    
    # get nullspace projected gradient
    i = n_sf;   j = n_sf - n_eq - ph.n_eq_off[1];

    # cannot use mul! to be allocation free here (variable Nullspace size) -> use old school style
    for k = 1:j
        ph.v_nem[k] = 0.0;
        for l = 1:i 
            ph.v_nem[k] += ph.dGdsf[l]*gm.N[l,k]
        end
    end

    for k = 1:i
        ph.grad[k] = 0.0;
        for l = 1:j 
            ph.grad[k] += ph.v_nem[l]*gm.N[k,l]
        end
    end
    

    return nothing

end

