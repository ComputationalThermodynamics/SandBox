# Amphibole data (Igneous database, Holland et al., 2018)
using StaticArrays;
using LinearAlgebra;
using BenchmarkTools


include("hb_pc_list.jl")

"""
    Structure to store solution phase data
"""
struct solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W,  
                    n_ox_t_n_em, n_sf_t_n_sf, n_eq_t_n_sf, n_em_t_n_W, n_em_t_n_sf,    # number of entries in matrixes
                    N, _C, _T, _I}
    ph          ::  String       # name of phase
    n_eq_off    ::  MVector{n_sf, _I}
    P           ::  _T
    T           ::  _T

    bulk        ::  MVector{n_ox, _T} 
    gb          ::  MVector{n_ox, _T}
    g0          ::  MVector{n_ox, _T}
    v_nem       ::  MVector{n_em, _T}
    v_nsf       ::  MVector{n_sf, _T}
    v_E         ::  MMatrix{n_ox, n_em, _T, n_ox_t_n_em}
    v_A         ::  MMatrix{n_sf, n_sf, _T, n_sf_t_n_sf}
    gamma       ::  MVector{n_em, _T}  

    emC         ::  MMatrix{n_ox, n_em, _T, n_ox_t_n_em}
    A           ::  MMatrix{n_eq, n_sf, _T, n_eq_t_n_sf}
    v           ::  MVector{n_ox,_T}
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

    ph          = "clino-amphibole"; 
    n_ox        = 11
    n_sf        = 17;
    n_eq        = 7;
    n_eq_off    = MVector{n_sf}(zeros(Int64,n_sf));
    n_em        = 11;
    n_xeos      = 10;
    n_W         = 55;   # perhaps this can be computed somehow from the variables above?

    P           = 5.0;                                  # pressure kbar
    T           = 650.0 + 273.15;                       # temperature C

    bulk        = SVector{n_ox}([0.50081,    0.08690,    0.11670,    0.12144,    0.07783,    0.00215,   0.02498,    0.01006,    0.00467,    0.00010,    0.05436]);

    # reference Gibbs energy of em, and Gamma at which the equilibrium was computed (obtained using MAGEMin)
    g0          = SVector{n_em}([-13012.62073, -13235.27114, -13472.30496, -12644.70794, -12762.02635, -10496.70590, -11477.04324, -11155.59746, -11828.15800, -13495.08535, -13063.17373]);
    gamma_eq    = SVector{n_ox}([-960.9655,    -1768.2476,   -788.4474,    -678.9683,    -355.2975,    -914.9708,    -839.9561,    -1008.3630,   -263.7269,    -1262.6087,   -368.4674]);
    gamma_em    = SVector{n_ox}([-960.965530,  -1747.660357, -795.996358,  -672.932965,  -348.636158,  -948.740786,  -851.768740,  -1014.787822, -245.521844,  -1115.578753, -368.238947]);
    gamma_pc    = SVector{n_ox}([-961.71227941038239350746, -1771.52191784699380150414, -783.52863672690375551610, -677.79113957713559557305, -353.68501121475600257327, -914.23627998836889219092, -832.35986650281336096668, -1012.53500402691236104147, -266.50224273504522898293, -1261.40041144356951008376, -368.23894675409849241987]);

    emC         = SMatrix{n_ox, n_em}(
                   [    +8.0  +0.0  +2.0  +5.0  +0.0  +0.0  +0.0  +0.0  +0.0  +0.0  +1.0;
                        +6.0  +2.0  +2.0  +3.0  +0.0  +0.0  +0.0  +0.0  +0.0  +0.0  +1.0;
                        +6.0  +1.5  +2.0  +4.0  +0.0  +0.0  +0.5  +0.0  +0.0  +0.0  +1.0;
                        +8.0  +1.0  +0.0  +3.0  +0.0  +0.0  +1.0  +0.0  +0.0  +0.0  +1.0;
                        +8.0  +0.0  +0.0  +7.0  +0.0  +0.0  +0.0  +0.0  +0.0  +0.0  +1.0;
                        +8.0  +0.0  +0.0  +0.0  +7.0  +0.0  +0.0  +0.0  +0.0  +0.0  +1.0;
                        +8.0  +0.0  +0.0  +3.0  +4.0  +0.0  +0.0  +0.0  +0.0  +0.0  +1.0;
                        +8.0  +0.0  +0.0  +2.0  +5.0  +0.0  +0.0  +0.0  +0.0  +0.0  +1.0;
                        +8.0  +0.0  +0.0  +3.0  +2.0  +0.0  +1.0  +0.0  +1.0  +0.0  +1.0;
                        +6.0  +1.5  +2.0  +4.0  +0.0  +0.5  +0.0  +0.0  +0.0  +0.0  +1.0;
                        +6.0  +1.0  +2.0  +3.0  +0.0  +0.0  +0.0  +2.0  +0.0  +0.0  +0.0]
                        );


    A           = SMatrix{n_eq, n_sf}(
                  [ 1 1 1 0 0    0    0  0  0 0    0    0    0 0 0 0 0;
                    0 0 0 1 1    0    0  0  0 0    0    0    0 0 0 0 0;
                    0 0 0 0 0    1    1  1  1 1    0    0    0 0 0 0 0;
                    0 0 0 0 0    0    0  0  0 0    1    1    1 1 0 0 0;
                    -1/4 0 0 0 0 -1/2 -1/2  0  0 0  1/2  1/2  1/2 0 1 0 0;
                    1/4 0 0 0 0  1/2  1/2  0  0 0 -1/2 -1/2 -1/2 0 0 1 0;
                    0 0 0 0 0   -1   -1 -1 -1 0    0    0    0 0 0 0 1]
                    );

    v            = SVector{n_em}([1,	1.5,	1.7,	0.80,	1,	1,	1,	1,	0.80,	1.7,	1.5]);
    W            = SVector{n_W}([20,	25,	65,	45,	75,	57,	63,	52,	30,	85,	-40, 25, 70, 80, 70, 72.50,	20,	-40, 35, 50, 90, 106.7,	94.80,	94.80, 40, 8, 15, 100, 113.5, 100,	111.2,	0,	54,	75,	33,	18,	23,	80,	87,	100,	12,	8,	91,	96,	65,	20,	80,	94,	95,	90,	94,	95,	50,	50,	35]);
    v_W          = MMatrix{n_em,n_W}(zeros(n_em,55));
    
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

    hb =  solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W, 
                         n_ox*n_em, n_sf*n_sf, n_eq*n_sf, n_em*n_W, n_em*n_sf,  
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
    get_gb!(hb);
    get_dpdsf!(hb);
    get_v_W!(hb);

    return hb;
end

"""
    calculate normalizing factor and needed partial derivatives
"""
function get_v_W!(ph::solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W}) where {n_ox, n_sf, n_eq, n_em, n_xeos, n_W}

    v_W = zeros(Float64, n_em, n_W)
    for i=1:n_em
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
    
    ph.p[1]  = ph.sf[11] - 2*ph.sf[16];
    ph.p[2]  = -ph.sf[2] - ph.sf[3] - ph.sf[10] + 2*ph.sf[16];
    ph.p[3]  = ph.sf[2];
    ph.p[4]  = -ph.sf[9] + ph.sf[14];
    ph.p[5]  = ph.sf[12];
    ph.p[6]  = 1 + 1/2 *ph.sf[2] + 1/2 *ph.sf[3] - ph.sf[4] - ph.sf[6] + ph.sf[11] + ph.sf[12] - 2*ph.sf[16];
    ph.p[7]  = ph.sf[4] - ph.sf[11] - ph.sf[12] - ph.sf[14];
    ph.p[8]  = -1/2 *ph.sf[2] - 1/2 *ph.sf[3] + ph.sf[6] - ph.sf[11] - ph.sf[12] + 2*ph.sf[16];
    ph.p[9]  = ph.sf[9];
    ph.p[10] = ph.sf[3];
    ph.p[11] = ph.sf[10];

    return nothing
end


"""
    Computes partial derivative of endmember fraction as function of site fraction
"""
function get_dpdsf!(ph::solution_phase{n_ox, n_sf, n_eq, n_em}) where {n_ox, n_sf, n_eq, n_em}
    ph.dpdsf[1,11]  = 1.0; ph.dpdsf[1,16] = -2.0;
    ph.dpdsf[2,2]   = -1.0; ph.dpdsf[2,3] = -1.0; ph.dpdsf[2,10] = -1.0; ph.dpdsf[2,16] = 2.0; 
    ph.dpdsf[3,2]   = 1.0;
    ph.dpdsf[4,9]   = -1.0;ph.dpdsf[4,14]  = 1.0;
    ph.dpdsf[5,12]  = 1.0;
    ph.dpdsf[6,2]   = 0.5; ph.dpdsf[6,3]  = 0.5; ph.dpdsf[6,4]  = -1.0; ph.dpdsf[6,6]  = -1.0; ph.dpdsf[6,11]  = 1.0; ph.dpdsf[6,12]  = 1.0; ph.dpdsf[6,16]  = -2.0;
    ph.dpdsf[7,4]   = 1.0; ph.dpdsf[7,11]  = -1.0; ph.dpdsf[7,12]  = -1.0; ph.dpdsf[7,14]  = -1.0; 
    ph.dpdsf[8,2]   = -0.5; ph.dpdsf[8,3]  = -0.5; ph.dpdsf[8,6]  = 1.0; ph.dpdsf[8,11]  = -1.0; ph.dpdsf[8,12]  = -1.0; ph.dpdsf[8,16]  = 2.0; 
    ph.dpdsf[9,9]   = 1.0;
    ph.dpdsf[10,3]  = 1.0;
    ph.dpdsf[11,10] = 1.0;

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
    
    # extract factors from matrix (and apply power only once)
    sf4_3 = ph.sf[4]*ph.sf[4]*ph.sf[4];
    sf5_3 = ph.sf[5]*ph.sf[5]*ph.sf[5];
    sf1 = ph.sf[1]
    sf2 = ph.sf[2]
    sf3 = ph.sf[3]
    
    sf6   = ph.sf[6]
    sf6_2 = sf6*sf6
    sf7_2 = ph.sf[7]*ph.sf[7]
    sf8   = ph.sf[8]
    sf8_2 = sf8*sf8
    sf9_2 = ph.sf[9]*ph.sf[9]
    sf10_2 = ph.sf[10]*ph.sf[10]
    sf11_2 = ph.sf[11]*ph.sf[11]
    sf12_2 = ph.sf[12]*ph.sf[12]
    
    sf13_2 = ph.sf[13]*ph.sf[13]
    sf14_2 = ph.sf[14]*ph.sf[14] 
    sf15   = ph.sf[15]
    sf15_05 = sqrt(sf15)
    sf16_05 = sqrt(ph.sf[16])
    sf17_2 = ph.sf[17]*ph.sf[17]

    ph.idm[1]  = sf1* sf4_3 *sf6_2 * sf11_2 *sf15 *sf17_2;
    ph.idm[2]  = 2.0*sf1* sf4_3* sf8_2 *sf11_2 *sf15_05 *sf16_05 *sf17_2;
    ph.idm[3]  = 8.0*sf2* sf4_3* sf6 *sf8 *sf11_2 *sf15_05* sf16_05* sf17_2;
    ph.idm[4]  = sf1* sf4_3 *sf8_2 *sf14_2 *sf15 *sf17_2;
    ph.idm[5]  = sf1* sf4_3 *sf6_2 *sf12_2 *sf15 *sf17_2;
    ph.idm[6]  = sf1* sf5_3 *sf7_2 *sf13_2 *sf15 *sf17_2;
    ph.idm[7]  = sf1* sf4_3 *sf7_2 *sf13_2 *sf15 *sf17_2;
    ph.idm[8]  = sf1* sf5_3 *sf6_2 *sf13_2 *sf15 *sf17_2;
    ph.idm[9]  = sf1* sf4_3 *sf9_2 *sf14_2 *sf15 *sf17_2;
    ph.idm[10] = 8.0*sf3* sf4_3* sf6 *sf8 *sf11_2 *sf15_05 *sf16_05 *sf17_2;
    ph.idm[11] = 2.0*sf1* sf4_3* sf10_2 *sf11_2 *sf15_05 *sf16_05 *sf10_2;

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

    ph.cv[1]  = (3*ph.sf[5] + 2*ph.sf[7] + 2*ph.sf[13])/(3*ph.sf[5] + 2*ph.sf[7] + 2*ph.sf[13] + 3*ph.sf[4] + 2*ph.sf[6] + 2*ph.sf[12]);
    ph.cv[2]  = ph.sf[8];
    ph.cv[3]  = ph.sf[14];
    ph.cv[4]  = ph.sf[3] + ph.sf[2];
    ph.cv[5]  = ph.sf[3]/(ph.sf[3] + ph.sf[2]);
    ph.cv[6]  = ph.sf[11];
    ph.cv[7]  = ph.sf[9];
    ph.cv[8]  = ph.sf[10];
    ph.cv[9]  = ph.cv[1] - ph.sf[5]/(ph.sf[5] + ph.sf[4]);
    ph.cv[10] = ph.cv[1] - ph.sf[7]/(ph.sf[7] + ph.sf[6]);

    return nothing;
end


"""
    convert compositional variables to site fractions
    the original compositional variables names and site-fraction are used in order to make the code more readable
"""
function get_ig(cv)

    x  =  cv[1];
    y  =  cv[2];
    z  =  cv[3];
    a  =  cv[4];
    k  =  cv[5];
    c  =  cv[6];
    f  =  cv[7];
    t  =  cv[8];
    Q1 =  cv[9];
    Q2 =  cv[10];
    
    xvA    = 1 - a;
    xNaA   = a + (-a) *k;
    xKA    = a* k;
    xMgM13 = 1 + Q1 - x;
    xFeM13 = -Q1 + x;
    xMgM2  = 1 - f + Q2 - t - x - y + (-f)* Q2 + (-Q2)* t + f* x + t* x + (-Q2)* y + x* y;
    xFeM2  = -Q2 + x + f* Q2 + Q2* t + (-f)* x + (-t)* x + Q2* y + (-x)* y;
    xAlM2  = y;
    xFe3M2 = f;
    xTiM2  = t;
    xCaM4  = c;
    xMgM4  = 1 - c - Q2 - x - z - (3/2)*Q1 + f*Q2 + Q2*t + c*x + Q2*y + x*z;
    xFeM4  = Q2 + x + 3/2* Q1 + (-f)* Q2 + (-Q2)* t + (-c)* x + (-Q2) *y + (-x) *z;
    xNaM4  = z;
    xSiT1  = 1 - 1/2* f - 1/2 *t - 1/2 *y + 1/2* z - 1/4* a;
    xAlT1  = 1/2* f + 1/2* t + 1/2 *y - 1/2* z + 1/4* a;
    xOHV   = 1 - t;

    #site fractions
    sf = [xvA,xNaA,xKA,xMgM13,xFeM13,xMgM2,xFeM2,xAlM2,xFe3M2,xTiM2,xCaM4,xMgM4,xFeM4,xNaM4,xSiT1,xAlT1,xOHV];

    return sf
end



"""
    Computes Gibbs energy and first derivative of the solution phase
"""
function compute_G!(ph,gv)

    get_p!(ph);
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
    # total time of routine: 1.505 μs -> 1.163 us
    get_p!(ph);         # 6.673 ns
    get_f!(ph,gv);      # 54.046 ns 
    get_Gex!(ph);       # 877.056 ns -> 759.440 ns
    get_idm!(ph);       # 212.879 ns  -> 102.421 ns -> 18.878 ns
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