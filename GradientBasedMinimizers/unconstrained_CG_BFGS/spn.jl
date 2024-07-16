# Spinel data (Igneous database, Holland et al., 2018)
using StaticArrays;
using LinearAlgebra;
using BenchmarkTools


include("spn_pc_list.jl")

"""
    Structure to store solution phase data
"""
struct solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W,  
                    n_ox_t_n_em, n_sf_t_n_sf, n_em_t_n_em, n_eq_t_n_sf, n_em_t_n_sf,    # number of entries in matrixes
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
    v           ::  MVector{1,_T}
    W           ::  MVector{n_W,_T}

    ig          ::  MVector{n_sf,_T}
    sf          ::  MVector{n_sf,_T}
    cv          ::  MVector{n_xeos,_T}
    p           ::  MVector{n_em,_T}
    f           ::  MVector{1,_T}
    df          ::  MVector{n_em,_T}
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

    ph          = "spinel"; 
    n_ox        = 11
    n_sf        = 10;
    n_eq        = 3;
    n_eq_off    = MVector{n_sf}(zeros(Int64,n_sf));
    n_em        = 8;
    n_xeos      = 7;
    n_W         = 28;   # perhaps this can be computed somehow from the variables above?

    bulk        = MVector{n_ox}([0.38494,    0.01776,    0.02824,    0.50566,    0.05886,    0.0001,   0.00250,    0.0010,    0.00096,    0.00109,    0.0]);
    bulk       .= bulk./(sum(bulk));

    # reference Gibbs energy of em, and Gamma at which the equilibrium was computed (obtained using MAGEMin)
    if test == 1
        P           = 12.0;                                  # pressure kbar
        T           = 1100.0 + 273.15;                       # temperature C    
        g0          = SVector{n_em}([-2515.94540,  -2500.25887,  -2217.27620, -2201.58966, -1452.03460, -1459.64806, -2033.47165, -2445.48343]);
        gamma_eq    = SVector{n_ox}([-1011.909631, -1829.092564, -819.264126, -695.467358, -412.948568, -971.890270, -876.544354, -1073.640927, -276.590707, -1380.299631, 0.0]);
        gamma_em    = SVector{n_ox}([-1016.399486, -1823.923434, -805.068209, -692.021970, -383.138855, -897.240517, -828.503925, -1062.477090, -310.231496, -1341.449679, 0.0]);
        gamma_pc    = SVector{n_ox}([-1006.885527, -1833.130473, -821.014905, -697.963191, -415.466554, -947.276954, -901.557511, -1058.244567, -236.176855, -1359.230803, 0.0]);
    elseif test == 2 #with solvi
        P           = 3.26;                                  # pressure kbar
        T           = 906.25 + 273.15;                       # temperature C    
        g0          = SVector{n_em}([-2491.74778,  -2474.94466,  -2185.44792, -2168.64480, -1408.67240, -1415.16928, -2004.02769, -2418.68444]);
        gamma_eq    = SVector{n_ox}([-1001.730935,-1818.611331,-812.972365,-689.113013,-396.911228,-966.511310,-882.719670,-1045.994137,-249.181839,-1332.815844,0.0]);
    end


    emC         = SMatrix{n_ox, n_em}(
                   [    +0.0  +1.0  +0.0  +1.0  +0.0  +0.0  +0.0  +0.0  +0.0  +0.0  +0.0
                        +0.0  +1.0  +0.0  +1.0  +0.0  +0.0  +0.0  +0.0  +0.0  +0.0  +0.0;
                        +0.0  +1.0  +0.0  +0.0  +1.0  +0.0  +0.0  +0.0  +0.0  +0.0  +0.0;
                        +0.0  +1.0  +0.0  +0.0  +1.0  +0.0  +0.0  +0.0  +0.0  +0.0  +0.0;
                        +0.0  +0.0  +0.0  +0.0  +3.0  +0.0  +0.0  +0.0  +1.0  +0.0  +0.0;
                        +0.0  +0.0  +0.0  +0.0  +3.0  +0.0  +0.0  +0.0  +1.0  +0.0  +0.0;
                        +0.0  +0.0  +0.0  +1.0  +0.0  +0.0  +0.0  +0.0  +0.0  +1.0  +0.0;
                        +0.0  +0.0  +0.0  +2.0  +0.0  +0.0  +0.0  +1.0  +0.0  +0.0  +0.0]
                        );


    A           = SMatrix{n_eq, n_sf}(
                    [       1.   1.  1. 1.  0.  0. 0. 0. 0. 0.;
                            0.5  0.5 0. 0.  2.  2. 1. 1. 1. 0.;
                           -0.5 -0.5 0. 0. -1. -1. 0. 0. 0. 1.]
                    );

    v            = SVector{1}([0.]);
    W            = SVector{n_W}([-8.2, 3.5, -13.0, 43.2, 49.1, -5.0, 22.5, 4.4, -6.0,36.8, 20.0, 14.0, 21.5, -8.2, 18.1, 49.0, -19.0, 35.1, -4.0, 7.6, -11.0, 9.0, 18.1, 11.9, 62.2, -6.4, 24.3, 60.0]);

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

    spn =  solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W, 
                         n_ox*n_em, n_sf*n_sf, n_em*n_em, n_eq*n_sf, n_em*n_sf,  
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

                            ig,
                            sf,
                            cv,
                            p,
                            f,
                            df,

                            Gex,
                            idm,
                            mu,
                            Graw,
                            G,
                            dpdsf,
                            dGdsf,
                            grad,
                            sf_off,
                            pcs     );

    # change of base with respect to the Gibbs hyperplane given
    # this function update g0 and the subsequent minimization is achieved with respect to the Gibbs hyperplane given by gamma
    get_gb!(spn);

    return spn;
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

    ph.p[1]          = -ph.sf[9] - 2*ph.sf[5]/3 + 2*ph.sf[1]/3 - 2*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - 4*ph.sf[10]/3 - (2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) + 1/3;
    ph.p[2]          = 2*ph.sf[5]/3 - 2*ph.sf[1]/3 - 4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - 2*ph.sf[10]/3 - 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) + 2/3;
    ph.p[3]          = ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - 2*ph.sf[8]/3 + 2*ph.sf[4]/3 - 2*ph.sf[6]/3 + 2*ph.sf[2]/3 + 2*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 2*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - (2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + (2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));
    ph.p[4]          = 2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 2*ph.sf[8]/3 - 2*ph.sf[4]/3 + 2*ph.sf[6]/3 - 2*ph.sf[2]/3 + 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));
    ph.p[5]          = -ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 2*ph.sf[8]/3 - 2*ph.sf[4]/3 - 2*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + (2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));
    ph.p[6]          = -2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - 2*ph.sf[8]/3 + 2*ph.sf[4]/3 - 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));
    ph.p[7]          = ph.sf[9];
    ph.p[8]          = 2*ph.sf[10];

    return nothing
end


"""
    Computes partial derivative of endmember fraction as function of site fraction
"""
function get_dpdsf!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}


    ph.dpdsf[1,1] = 2*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + (2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 2/3;      ph.dpdsf[1,2] = 2*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 2*ph.sf[10]/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) + (2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 1/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      ph.dpdsf[1,3] = 0;      ph.dpdsf[1,4] = 0;      ph.dpdsf[1,5] = 4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 2/3;      ph.dpdsf[1,6] = 4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 4*ph.sf[10]/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) + 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 2/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      ph.dpdsf[1,7] = 0;      ph.dpdsf[1,8] = 0;      ph.dpdsf[1,9] = -1;      ph.dpdsf[1,10] = -2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - 4/3;      
    ph.dpdsf[2,1] = 4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 2/3;      ph.dpdsf[2,2] = 4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 4*ph.sf[10]/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) + 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 2/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      ph.dpdsf[2,3] = 0;      ph.dpdsf[2,4] = 0;      ph.dpdsf[2,5] = 8*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 4*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 2/3;      ph.dpdsf[2,6] = 8*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 8*ph.sf[10]/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) + 4*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 4/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      ph.dpdsf[2,7] = 0;      ph.dpdsf[2,8] = 0;      ph.dpdsf[2,9] = 0;      ph.dpdsf[2,10] = -4*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - 2/3;      
    ph.dpdsf[3,1] = -2*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - (2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2);      ph.dpdsf[3,2] = -2*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 2*ph.sf[10]/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - (2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 2/3 + 1/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      ph.dpdsf[3,3] = -ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + (2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2);      ph.dpdsf[3,4] = -ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + ph.sf[9]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - 2*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2*ph.sf[10]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + (2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2/3 - 1/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[3,5] = -4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2);      ph.dpdsf[3,6] = -4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 4*ph.sf[10]/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 2/3 + 2/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      ph.dpdsf[3,7] = -2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2);      ph.dpdsf[3,8] = -2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2*ph.sf[9]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 4*ph.sf[10]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2/3 - 2/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[3,9] = (2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[3,10] = 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      
    ph.dpdsf[4,1] = -4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2);      ph.dpdsf[4,2] = -4*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 4*ph.sf[10]/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - 2*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 2/3 + 2/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      ph.dpdsf[4,3] = -2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2);      ph.dpdsf[4,4] = -2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2*ph.sf[9]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 4*ph.sf[10]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2/3 - 2/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[4,5] = -8*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) - 4*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2);      ph.dpdsf[4,6] = -8*ph.sf[10]*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 8*ph.sf[10]/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])) - 4*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1])^2) + 2/3 + 4/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      ph.dpdsf[4,7] = -4*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 8*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 4*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2);      ph.dpdsf[4,8] = -4*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 4*ph.sf[9]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - 8*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 8*ph.sf[10]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 4*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2/3 - 4/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[4,9] = 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[4,10] = 4*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 4*(2*ph.sf[6] + ph.sf[2])/(3*(2*ph.sf[6] + ph.sf[2] + 2*ph.sf[5] + ph.sf[1]));      
    ph.dpdsf[5,1] = 0;      ph.dpdsf[5,2] = 0;      ph.dpdsf[5,3] = ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - (2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2);      ph.dpdsf[5,4] = ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - ph.sf[9]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 2*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2*ph.sf[10]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - (2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2/3 + 1/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[5,5] = 0;      ph.dpdsf[5,6] = 0;      ph.dpdsf[5,7] = 2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2);      ph.dpdsf[5,8] = 2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2*ph.sf[9]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 4*ph.sf[10]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2/3 + 2/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[5,9] = -(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[5,10] = -2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      
    ph.dpdsf[6,1] = 0;      ph.dpdsf[6,2] = 0;      ph.dpdsf[6,3] = 2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2);      ph.dpdsf[6,4] = 2*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2*ph.sf[9]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 4*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 4*ph.sf[10]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - 2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 2/3 + 2/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[6,5] = 0;      ph.dpdsf[6,6] = 0;      ph.dpdsf[6,7] = 4*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) + 8*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 4*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2);      ph.dpdsf[6,8] = 4*ph.sf[9]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 4*ph.sf[9]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) + 8*ph.sf[10]*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 8*ph.sf[10]/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])) - 4*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4])^2) - 2/3 + 4/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[6,9] = -2*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      ph.dpdsf[6,10] = -4*(2*ph.sf[8] + ph.sf[4])/(3*(2*ph.sf[7] + ph.sf[3] + 2*ph.sf[8] + ph.sf[4]));      
    ph.dpdsf[7,1] = 0;      ph.dpdsf[7,2] = 0;      ph.dpdsf[7,3] = 0;      ph.dpdsf[7,4] = 0;      ph.dpdsf[7,5] = 0;      ph.dpdsf[7,6] = 0;      ph.dpdsf[7,7] = 0;      ph.dpdsf[7,8] = 0;      ph.dpdsf[7,9] = 1;      ph.dpdsf[7,10] = 0;      
    ph.dpdsf[8,1] = 0;      ph.dpdsf[8,2] = 0;      ph.dpdsf[8,3] = 0;      ph.dpdsf[8,4] = 0;      ph.dpdsf[8,5] = 0;      ph.dpdsf[8,6] = 0;      ph.dpdsf[8,7] = 0;      ph.dpdsf[8,8] = 0;      ph.dpdsf[8,9] = 0;      ph.dpdsf[8,10] = 2;  
 
    return nothing
end

"""
    Calculate G excess
"""
function get_Gex!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}

    a = ph.v_E;
    b = ph.p;
    c = ph.W;

    @inbounds for i=1:n_em
        it      = 1;
        Gex_v = 0.0;
        @inbounds for j=1:n_em-1
            tmp = (a[i,j]-b[j]);
            for k=j+1:n_em
                Gex_v -= tmp*(a[i,k]-b[k])*c[it];
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

    sf5 = sqrt(ph.sf[5]);
    sf6 = sqrt(ph.sf[6]);
    sf7 = sqrt(ph.sf[7]);

    ph.idm[1]          = ph.sf[7]*ph.sf[1];
    ph.idm[2]          = 2*sf7*ph.sf[3]*sf5;
    ph.idm[3]          = ph.sf[7]*ph.sf[2];
    ph.idm[4]          = 2*sf7*ph.sf[3]*sf6;
    ph.idm[5]          = ph.sf[8]*ph.sf[2];
    ph.idm[6]          = 2*sqrt(ph.sf[8])*ph.sf[4]*sf6;
    ph.idm[7]          = ph.sf[9]*ph.sf[1];
    ph.idm[8]          = 2*sf5*ph.sf[1]*sqrt(ph.sf[10]);

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
    #mul!(ph.v_nem,ph.emC,ph.gamma);
    ph.v_nem .= ph.emC*ph.gamma
    ph.gb .= ph.g0 .- ph.v_nem;
    return nothing
end



function get_cv!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}

    ph.cv[1] = (2 *ph.sf[6] + ph.sf[2])/(2 *ph.sf[6] + ph.sf[2] + 2 *ph.sf[5] + ph.sf[1])
    ph.cv[2] = (2 *ph.sf[8] + ph.sf[4])/(2 *ph.sf[7] + ph.sf[3] + 2 *ph.sf[8] + ph.sf[4])
    ph.cv[3] = ph.sf[9]
    ph.cv[4] = 2 *ph.sf[10]
    ph.cv[5] = -ph.sf[5] + ph.sf[1]
    ph.cv[6] = -ph.sf[6] + ph.sf[2]
    ph.cv[7] = -ph.sf[8] + ph.sf[4]

    return nothing;
end



"""
    convert compositional variables to site fractions
    the original compositional variables names and site-fraction are used in order to make the code more readable
"""
function get_ig(cv)

    x  =  cv[1];
    y  =  cv[2];
    c  =  cv[3];
    t  =  cv[4];
    Q1 =  cv[5];
    Q2 =  cv[6];
    Q3 =  cv[7];

    xMgT = 1/3 + 1/3*t - 1/3*x + 2/3*Q1 + (-1/3*t)*x;
    xFeT = 1/3*x + 2/3*Q2 + 1/3*t*x;
    xAlT = 2/3 - 1/3*t - 2/3*Q1 - 2/3*Q2 - 2/3*Q3 - 2/3*y + 2/3*c*y + 2/3*t*y;
    xFe3T = 2/3*Q3 + 2/3*y + (-2/3*c)*y + (-2/3*t)*y;
    xMgM = 1/3 - 1/3*Q1 + 1/3*t - 1/3*x + (-1/3*t)*x;
    xFeM = -1/3*Q2 + 1/3*x + 1/3*t*x;
    xAlM = 2/3 + 1/3*Q1 + 1/3*Q2 + 1/3*Q3 - c - 2/3*y - 5/6*t + 2/3*c*y + 2/3*t*y;
    xFe3M = -1/3*Q3 + 2/3*y + (-2/3*c)*y + (-2/3*t)*y;
    xCrM = c;
    xTiM = 1/2*t;

    #site fractions
    sf = [xMgT,xFeT,xAlT,xFe3T,xMgM,xFeM,xAlM,xFe3M,xCrM,xTiM];

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
    # total time of routine: 1.505 μs
    get_p!(ph);         # 6.673 ns
    get_dpdsf!(ph);
    get_f!(ph,gv);      # 54.046 ns 
    get_Gex!(ph);       # 877.056 ns -> 759.440 ns
    get_idm!(ph);       # 212.879 ns  -> 102.421 ns 
    get_mu!(ph,gv);     # 72.523 ns

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

