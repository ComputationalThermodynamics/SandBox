# Gradient based method: data and routines
struct gradient_method{n_sf, n_sf_t_n_sf, _T,_I}

    G_track  ::  MVector{2048,_T}
    G0_out   ::  MVector{1,_T}
    G0       ::  MVector{1,_T}

    ite      ::  MVector{1,_I}
    i        ::  MVector{1,_I}
    o        ::  MVector{1,_I}

    x0       ::  MVector{n_sf,_T}
    x        ::  MVector{n_sf,_T}

    Ac       ::  MVector{n_sf,_T}
    Bc       ::  MVector{n_sf,_T}
    Cc       ::  MVector{n_sf,_T}

    xS       :: MMatrix{64, n_sf,_T} 

    grad0    ::  MVector{n_sf,_T}
    grad1    ::  MVector{n_sf,_T}

    sk       ::  MVector{n_sf,_T}
    yk       ::  MVector{n_sf,_T}
    t1       ::  MMatrix{n_sf,n_sf,_T, n_sf_t_n_sf} 
    t2       ::  MMatrix{n_sf,n_sf,_T, n_sf_t_n_sf} 
    v_T      ::  MMatrix{n_sf,n_sf,_T, n_sf_t_n_sf} 
    v_T2     ::  MMatrix{n_sf,n_sf,_T, n_sf_t_n_sf} 
    v_t      ::  MVector{n_sf,_T}

    BkI      ::  MMatrix{n_sf,n_sf,_T, n_sf_t_n_sf} 
    N        ::  MMatrix{n_sf,n_sf,_T, n_sf_t_n_sf} 
    pk       ::  MVector{n_sf,_T}
    pkur     ::  MVector{n_sf,_T}

    ur       ::  MVector{1,_T}
    urm      ::  MVector{1,_T}
    alpha    ::  MVector{1,_T}

    dxi      ::  MVector{1,_T}
    dxo      ::  MVector{1,_T}
    beta_cg  ::  MVector{1,_T}
    beta     ::  _T
    eta      ::  _T
    tol_dx   ::  _T
    eps      ::  _T
    bnd      ::  _T
    max_ite  ::  _I
end

"""
    function to allocate the memory necessary to conducte gradient based optimization methods
"""
function init_gradient_methods(ph::solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W}) where {n_ox, n_sf, n_eq, n_em, n_xeos, n_W}
    G_track  = MVector{2048}(zeros(2048));
    G0_out    = MVector{1}(zeros(1));
    G0        = MVector{1}(zeros(1));

    ite       = MVector{1}(zeros(Int64,1));   
    i         = MVector{1}(zeros(Int64,1));   
    o         = MVector{1}(zeros(Int64,1));   

    x0        = MVector{n_sf}(zeros(n_sf));
    x         = MVector{n_sf}(zeros(n_sf));

    Ac        = MVector{n_sf}(zeros(n_sf));
    Bc        = MVector{n_sf}(zeros(n_sf));
    Cc        = MVector{n_sf}(zeros(n_sf));

    xS        = MMatrix{64, n_sf}(zeros(64,n_sf));

    grad0     = MVector{n_sf}(zeros(n_sf));
    grad1     = MVector{n_sf}(zeros(n_sf));

    sk        = MVector{n_sf}(zeros(n_sf));
    yk        = MVector{n_sf}(zeros(n_sf));
    t1        = MMatrix{n_sf, n_sf}(zeros(n_sf,n_sf));
    t2        = MMatrix{n_sf, n_sf}(zeros(n_sf,n_sf));
    v_T       = MMatrix{n_sf, n_sf}(zeros(n_sf,n_sf));
    v_T2      = MMatrix{n_sf, n_sf}(zeros(n_sf,n_sf));
    v_t       = MVector{n_sf}(zeros(n_sf));

    BkI       = MMatrix{n_sf, n_sf}(Matrix(1.0I, n_sf,n_sf))
    Ni        = MMatrix{n_sf, n_sf}(zeros(n_sf,n_sf));
    
    pk        = MVector{n_sf}(zeros(n_sf));
    pkur      = MVector{n_sf}(zeros(n_sf));
    alpha     = MVector{1}(zeros(1));
    ur        = MVector{1}(zeros(1));
    urm       = MVector{1}(zeros(1));
    dxi       = MVector{1}(zeros(1));
    dxo       = MVector{1}(zeros(1));
    beta_cg   = MVector{1}(zeros(1));

    beta      = 0.25;
    eta       = 1e-3;
    tol_dx    = 1e-7;
    eps       = 1e-9;
    bnd       = 1e-8;

    max_ite   = 256;

    gm =  gradient_method{n_sf, n_sf*n_sf, Float64, Int64}(  
                                G_track,    
                                G0_out,
                                G0,

                                ite,
                                i,
                                o,
        
                                x0,
                                x,

                                Ac,
                                Bc,
                                Cc,

                                xS,

                                grad0,
                                grad1,

                                sk,
                                yk,
                                t1,
                                t2,
                                v_T,
                                v_T2,
                                v_t,

                                BkI,
                                Ni,
                                pk,
                                pkur,

                                ur,
                                urm,
                                alpha,
                                dxi,
                                dxo,

                                beta_cg,
                                beta,
                                eta,
                                tol_dx,
                                eps,
                                bnd,
                                
                                max_ite     )

    return gm
end


"""
    function to update pk, the descent direction
"""
function   update_pk_BFGS!(gv,ph,gm)

    compute_G_dG!(gm,ph,gv);
    gm.G0[1]     = ph.G[1];
    gm.grad0    .= ph.grad;

    gm.pk       .= gm.BkI*ph.grad;

    # get nullspace projected descent direction. This part is critical to ensure that the descent direction respect at all time the equality constraints
    i = length(ph.sf);   j = i - size(ph.A,1) - ph.n_eq_off[1];

    # cannot use mul! to be allocation free here (variable Nullspace size) -> use old school style
    for k = 1:j
        ph.v_nem[k] = 0.0;
        for l = 1:i 
            ph.v_nem[k] += gm.pk[l]*gm.N[l,k]
        end
    end

    for k = 1:i
        gm.pk[k] = 0.0;
        for l = 1:j 
            gm.pk[k] += ph.v_nem[l]*gm.N[k,l]
        end
    end

    return nothing

end

"""
    function to update pk, the descent direction
"""
function   update_pk_BFGS_CG!(gv,ph,gm)

    compute_G_dG!(gm,ph,gv);
    gm.G0[1]     = ph.G[1];
    gm.grad0    .= ph.grad;

    if (gm.i == 1)
        gm.pk .= gm.BkI*ph.grad;
    else
        gm.pk .= gm.BkI*ph.grad .+ gm.eta*gm.beta_cg[1]*gm.pk;
    end
    # get nullspace projected descent direction. This part is critical to ensure that the descent direction respect at all time the equality constraints
    i = length(ph.sf);   j = i - size(ph.A,1) - ph.n_eq_off[1];

    # cannot use mul! to be allocation free here (variable Nullspace size) -> use old school style
    for k = 1:j
        ph.v_nem[k] = 0.0;
        for l = 1:i 
            ph.v_nem[k] += gm.pk[l]*gm.N[l,k]
        end
    end

    for k = 1:i
        gm.pk[k] = 0.0;
        for l = 1:j 
            gm.pk[k] += ph.v_nem[l]*gm.N[k,l]
        end
    end


    return nothing

end

"""
    function to retrieve the maximum allowed step before violating bounds
"""
function get_ur_dist!(ph,gm)

    dl      = 1.0;      dl0     = 1.0;
                        dlm0    = 1e10;
     
    for i=1:lastindex(ph.sf)
        if ph.sf[i] < gm.eps
            dl  = (gm.x[i] - gm.eps*5.0)/abs(gm.pk[i]);

            if dl < dl0
                dl0 = dl;
            end

        end
        if gm.pk[i] > 0.0 #if descent direction decreases site fraction
            dl  = (gm.x[i] - gm.eps*5.0)/gm.pk[i];
            if dl < dlm0
                dlm0 = dl;
            end  
        end

    end

    gm.ur  .= 1.0/dl0;
    gm.urm .= 1.0/dlm0;

    return nothing
end


"""
    Backtracking line search using the Armijo condition
"""
function backTrackingLineSearch!(gv,ph,gm)

    G0          =  ph.G[1]
    # gm.alpha[1] =  gm.dxi[1]/100.0;
    gm.alpha[1] = 1e-4
    t           =  1.0;   

    gm.pkur    .= .-gm.pk./gm.ur[1];
    ph.sf      .=  gm.x .+ t.*gm.pkur;

    stp         = max(gm.pkur' * ph.grad, -1.0)

    compute_G!(ph,gv);
    lhs         =  ph.G[1];
    rhs         =  G0 + gm.alpha[1]*t*stp;

    while lhs > rhs
        t       =  gm.beta*t;
        ph.sf  .=  gm.x .+ t.*gm.pkur;
        compute_G!(ph,gv);
        lhs     =  ph.G[1];
        rhs     =  G0 + gm.alpha[1]*t*stp;
    end

    gm.x       .=  ph.sf;
    return nothing
end

"""
    Backtracking line search using the wolf condition
"""
function backTrackingLineSearchW!(gv,ph,gm)

    G0          =  ph.G[1]
    # gm.alpha[1] =  gm.dxi[1]/10.0;
    gm.alpha[1] =  0.001;
    t           =  1.0;   
    gm.pkur    .= .-gm.pk./gm.ur[1];
    ph.sf      .=  gm.x .+ t.*gm.pkur;

    stp         =  gm.pkur'*ph.grad;
    if stp < -1.0
        stp = -1.0
    end

    compute_G_dG!(gm,ph,gv);
    lhs         =  ph.G[1];
    rhs         =  G0 + gm.alpha[1]*t*stp;

    lhs2        = 0.0;
    rhs2        = 0.5*ph.grad'*gm.pkur;

    while lhs > rhs && lhs2 < rhs2
        t       =  gm.beta*t;
        ph.sf  .=  gm.x .+ t.*gm.pkur;
        compute_G_dG!(gm,ph,gv);
        lhs     =  ph.G[1];
        rhs     =  G0 + gm.alpha[1]*t*stp;
        lhs2    =  ph.grad'*gm.pkur;
    end

    gm.x       .=  ph.sf;

    return nothing
end

"""
    update pseudo-hessian inverse using Sherman-Morrison formulae 
"""
function update_BkI!(gv,ph,gm)

    G0          = gm.G0[1]; 

    ph.sf      .=  gm.x;
    compute_G_dG!(gm,ph,gv);
    G1          = ph.G[1];
    gm.grad1   .= ph.grad;

    gm.sk      .= gm.x .- gm.x0;     
    gm.yk      .= (gm.grad1 .- gm.grad0);

    if norm(gm.sk) > 0.0

        gm.v_T      .= gm.sk*gm.sk';
        gm.v_t      .= gm.BkI*gm.yk;
        tmp         = gm.sk'*gm.yk;
        gm.v_T    ./= (tmp^2);
        gm.t1      .= gm.v_T .* (tmp + gm.yk'*gm.v_t);

        gm.v_T     .= gm.v_t*gm.sk';
        gm.t2      .= gm.sk*gm.yk';
        gm.v_T2    .= gm.t2*gm.BkI;
        gm.t2      .= (gm.v_T .+ gm.v_T2)./tmp;

        gm.BkI    .= gm.BkI .+ (gm.t1 .- gm.t2);
    end
    gm.dxi    .= abs((G1-G0)/G0);

    return nothing
end


"""
    This function checks for the inactive site fractions ( bnd < value < tol ) and add them to the computation of the nullspace
    The function is called at the beginning of every outter loop to update the gradient based direction
"""
function update_Nullspace!(gm,ph::solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W}) where {n_ox, n_sf, n_eq, n_em, n_xeos, n_W}

    ph.n_eq_off[1]  = 0;

    ph.v_A[1:size(ph.A)[1],1:size(ph.A)[2]] = ph.A;

    for i=1:lastindex(ph.sf)
        if ph.sf_off[i] == 1
            ph.n_eq_off[1] += 1;
            ph.v_A[size(ph.A)[1]+ph.n_eq_off[1],1:size(ph.A)[2]] = ph.v_E[i,:]';
        end
    end
    
    #get the right size of the nullspace
    i = n_sf;   j = n_sf - n_eq - ph.n_eq_off[1];

    gm.N[1:i,1:j] = nullspace(ph.v_A[1:size(ph.A)[1]+ph.n_eq_off[1], 1:i]);  

    return nothing
end


"""
    reset BkI
"""
function reset_BkI!(gm,ph::solution_phase{n_ox, n_sf, n_eq, n_em, n_xeos, n_W}) where {n_ox, n_sf, n_eq, n_em, n_xeos, n_W}

    for i=1:n_sf
        for j=1:n_sf
            gm.BkI[i,j] = 0.0;
        end
        gm.BkI[i,i] = 1.0;
    end

    return nothing
end

"""
    Minimizer function: BFGS method
"""
function null_min_BFGS!(gv,ph,gm)

    ph.sf          .= ph.ig;
    gm.x           .= ph.sf;
    gm.x0          .= ph.sf;

    gm.G0_out[1]    = 1.0;
    gm.dxo[1]       = 1.0;             # initialize outter loop tolerance

    omax            = 128;
    imax            = 512;

    gm.o[1] = 1;     gm.ite[1] = 1; 
    while gm.dxo[1] > gm.tol_dx && gm.o[1] < omax
        gm.dxi[1]        = 1.0;             # initialize inner loop tolerance

        reset_BkI!(gm,ph); 

        gm.i[1] = 1;
        while gm.dxi[1] > gm.tol_dx && gm.i[1] < imax

            gm.x0       .= gm.x;            #
            ph.sf       .= gm.x;

            update_pk_BFGS!(gv,ph,gm);   
            ph.sf       .= gm.x .- gm.pk;

            get_ur_dist!(ph,gm);

            backTrackingLineSearch!(gv,ph,gm)  
                
            update_BkI!(gv,ph,gm);  

            gm.G_track[gm.ite[1]] = gm.dxi[1];

            gm.i[1]     += 1;
            gm.ite[1]   += 1;
        end

        gm.dxo[1]       = abs((ph.G[1]-gm.G0_out[1])/gm.G0_out[1]);
        gm.G0_out[1]    = ph.G[1];
        gm.o[1]        += 1;
    end

    return nothing
end

"""
    Minimizer function: BFGS_s method
    same method than BFGS, we just save the minimization path here (we are not interested in performances)
"""
function null_min_BFGS_s!(gv,ph,gm,trackSF) #

    ph.sf          .= ph.ig;
    gm.x           .= ph.sf;
    gm.x0          .= ph.sf;

    gm.G0_out[1]    = 1.0;
    gm.dxo[1]       = 1.0;             # initialize outter loop tolerance

    omax            = 128;
    imax            = 512;

    gm.o[1] = 1;     gm.ite[1] = 1; 
    while gm.dxo[1] > gm.tol_dx && gm.o[1] < omax
        gm.dxi[1]        = 1.0;             # initialize inner loop tolerance

        reset_BkI!(gm,ph); 

        gm.i[1] = 1;
        while gm.dxi[1] > gm.tol_dx && gm.i[1] < imax

            push!(trackSF,Array(gm.x));

            gm.x0       .= gm.x;            #
            ph.sf       .= gm.x;

            update_pk_BFGS!(gv,ph,gm);   
            ph.sf       .= gm.x .- gm.pk;

            get_ur_dist!(ph,gm);

            backTrackingLineSearch!(gv,ph,gm)  
                
            update_BkI!(gv,ph,gm);  

            gm.G_track[gm.ite[1]] = gm.dxi[1];

            gm.i[1]     += 1;
            gm.ite[1]   += 1;
        end

        gm.dxo[1]       = abs((ph.G[1]-gm.G0_out[1])/gm.G0_out[1]);
        gm.G0_out[1]    = ph.G[1];
        gm.o[1]        += 1;
    end

    return trackSF
end

"""
    Minimizer function: Conjugate gradient method
"""
function null_min_CG!(gv,ph,gm)
    gm.x           .= ph.ig;

    gm.G0_out[1]    = 1.0;
    gm.dxo[1]       = 1.0;                  # initialize outter loop tolerance
    G0              = 1.0; 

    omax            = 1024;
    imax            = 1024;

    gm.o[1] = 1;     gm.ite[1] = 1; 
    while gm.dxo[1] > gm.tol_dx && gm.o[1] < omax
        gm.dxi[1]        = 1.0;             # initialize inner loop tolerance

        # initialize gradient and descent direction
        ph.sf          .= gm.x;
        compute_G_dG!(gm,ph,gv);

        gm.pk          .= ph.grad;
        gm.grad0       .= ph.grad;

        gm.i[1] = 1;
        while gm.dxi[1] > gm.tol_dx && gm.i[1] < imax
            ph.sf       .= gm.x .- gm.pk;

            get_ur_dist!(ph,gm);

            backTrackingLineSearch!(gv,ph,gm)  
            
            ph.sf       .= gm.x;
            compute_G_dG!(gm,ph,gv);
            G1           = ph.G[1]; 
            gm.grad1    .= ph.grad;

            gm.sk       .= gm.pk.*gm.pk;
            n            = sum(gm.sk);

            gm.v_t      .= gm.grad1 - gm.grad0;
            gm.yk       .= gm.v_t + gm.pk;


            beta         = (gm.grad1'*gm.yk)./sqrt(n);                  # Beta from Rivaie et al., 2012 AMC

            # gm.sk       .= (gm.yk .+ gm.pk)
            # beta         = (gm.grad1'*gm.sk)./sqrt(n);                # Beta from Rivaie et al., 2015 AMC
            
            theta        = (gm.grad1'*gm.pk)./sqrt(n);                  # thetafrom Liu et al., 2018

            gm.pk       .= (gm.grad1 .+ beta*gm.pk .- theta*gm.v_t);


            # get nullspace projected descent direction. This part is critical to ensure that the descent direction respect at all time the equality constraints
            i = length(ph.sf);   j = i - size(ph.A,1) - ph.n_eq_off[1];

            # cannot use mul! to be allocation free here (variable Nullspace size) -> use old school style
            for k = 1:j
                ph.v_nem[k] = 0.0;
                for l = 1:i 
                    ph.v_nem[k] += gm.pk[l]*gm.N[l,k]
                end
            end

            for k = 1:i
                gm.pk[k] = 0.0;
                for l = 1:j 
                    gm.pk[k] += ph.v_nem[l]*gm.N[k,l]
                end
            end


            gm.grad0    .= gm.grad1;
            gm.dxi      .= abs((G1-G0)/G0);
            G0           = ph.G[1]; 

            gm.G_track[gm.ite[1]] = G1;

            gm.i[1]     += 1;
            gm.ite[1]   += 1;
        end

        gm.dxo[1]       = abs((G0-gm.G0_out[1])/gm.G0_out[1]);
        gm.G0_out[1]    = ph.G[1];
        gm.o[1]        += 1;
    end

    return nothing
end