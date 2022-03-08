%%
% Simple matlab script to present the minimization approach employed in MAGEMin (Riel et al., 2022)
% Plagioclase model after Holland et al., 2021
% Linear programming approach described in de Capitani & Brown (1987)

clear; close all;
set(0,'defaultfigurecolor',[1 1 1]);

%% control panel
%Run with LP, 0 or PGE, 1
PGE         = 1;
max_ite     = 24;

%fmincon solver options
options     = optimoptions('fmincon','Algorithm','sqp-legacy','OptimalityTolerance',1e-10,'Display', 'off');    % run interior-point algorithm

%% System definition
Gamma_sol = [-970.718847000000,-1770.76203800000,-761.944990621532,-837.091180964866,-767.348396813270];

% pressure and temperature conditions
P   = 3;                                                                    %kbar
T   = 600 + 273.15;                                                         %Kelvin
R   = 0.0083144;                                                            %gas constant

oxi = {'SiO2' 'Al2O3' 'CaO' 'K2O' 'Na2O'};                                  %NCKAS system
apo = [3 5 2 3 3];                                                          %number of atoms per oxide
b   = [0.7069    0.1663    0.0456    0.0445    0.0367];                     %bulk rock composition

%% Phases definition
% Quartz and Sillimanite definition 
gbqtz    = -970.718847;  natom_qtz = 3; facqtz   = sum(b.*apo)/natom_qtz;
gbsil    = -2741.480885; natom_sil = 8; facsil   = sum(b.*apo)/natom_sil;

Cqtz     = [+1        +0        +0         +0        +0        ];           %qtz
Csil     = [+1        +1        +0         +0        +0        ];           %sill

% pl4T definition 
n_em     = 3;                                                               % number of feldspar endmembers 
n_C      = 5;                                                               %number of oxides (system components)
n_pp     = 10;                                                              %number of pure phase to start with (includes 5 fakes oxide phases)

%endmember molar composition 	                                                 
C(1,1:5) = [+3.000000 +0.500000 +0.000000  +0.000000 +0.500000 ];           %ab
C(2,1:5) = [+2.000000 +1.000000 +1.000000  +0.000000 +0.000000 ];           %an
C(3,1:5) = [+3.000000 +0.500000 +0.000000  +0.500000 +0.000000 ];           %san

%normalisation facpl
natom_pl = sum(C(1,:).*apo);
facpl   = sum(b.*apo)/natom_pl;

% Gibbs energy of reference
gb(1)    = -4177.704449;
gb(2)    = -4470.047236;
gb(3)    = -4215.127702;

% plagioclase margules
W(1)     = 14.6 - 0.00935*T - 0.04*P;
W(2)     = 24.1 - 0.00957*T + 0.338*P;
W(3)     = 48.5 - 0.13*P;
v(1)     = 0.674;
v(2)     = 0.55;
v(3)     = 1.0;


%% calculate equilibrium
%allocate memory
Gamma       = zeros(1,n_C);                                                 %initialize Gamma
n_vec       = b;                                                            %initial fraction vector = bulk rock composition
idA         = zeros(1,n_C);
n_swap      = 0;

%generate list of discretized solution model points (pseudocompounds)
[G_,C_,x_]  = gen_PC_fct(Gamma,C,gb,v,W,T,R,facpl, gbsil,Csil,facsil, gbqtz,Cqtz,facqtz);
A(:,:)      = C_(1:n_C,1:n_C);                                              %initial composition = identity matrix
g_A(:)      = G_(1:n_C)';                                                   %initial gibbs energy = penalty

color(1:size(x_,1))     = 0.0;                                              %make the color of the minimized solution model depending on the iteration stage
xi_(1:size(x_,1),n_em)  = 0.0;
A1                      = inv(A);                                           %A is identity so A1 = A here, it is inversed for algorithm clarity
min_path                = [];

%% THERIAK-like minimization approach (linear programming)
if PGE == 0
    for ite=1:max_ite
        for lp=1:4
            for i=size(C,2):size(C_,1)
                B1      = A1*C_(i,:)';
                dG      = G_(i) - sum(B1.*g_A');

                ph2swp  = -1;
                if dG < 0.0
                    min_F  =  1e6;
                    for j=1:size(C,2)
                        F = n_vec(j)/B1(j);
                        if F < min_F && F > 0.0
                            min_F  = F;
                            ph2swp = j;
                        end
                    end
                end

                if ph2swp ~= -1
                    swp     = 1;
                    n_swap  = n_swap + 1;
                    idA(ph2swp) = i;
                    g_A(ph2swp)     = G_(i);

                    for j=1:size(C,2)
                        A(:,ph2swp) = C_(i,:)';
                    end

                    A0      = A;
                    A1      = inv(A0);   
                    n_vec   = A\b';
                end

            end  
        end
        Gamma           = (A'\g_A')';
        
        %calculate the residual with respect to the reference Gamma
        Residual(ite)   = norm(Gamma_sol-Gamma);

        rot = rotate_ghpp(C,Gamma);                                         % apply change of base to the pure endmembers
        xp = [];
        for i=1:n_C
            if idA(i) > n_pp                                                % if the pseudocompound does not belong to a fake phase or a pure endmember
                sol     = fmincon(@(x) calc_G_fct(rot,x,gb,v,W,T,R,facpl),[x_(idA(i),1) x_(idA(i),2)],[],[],[],[],[1e-8 1e-8],[1-1e-8 1-1e-8],[],options);
                [G,phC] = calc_GC_fct(zeros(1,size(C,2)),[sol(1) sol(2)],C,gb,v,W,T,R,facpl);
       
                G_      = [G_, G];
                C_      = [C_; phC];
                x_      = [x_; [sol(1) sol(2)]];
                color   = [color, min(0.0+(0.9/max_ite*3)*ite,1)];
                xp      = [xp, [sol(1) sol(2)]];
                [G,phC,xi]      = calc_GC_fct(rot,[sol(1) sol(2)],C,gb,v,W,T,R,facpl);
                xi_(idA(j),:)   = xi(1,:);
            end
        end
        min_path = [min_path; xp];
        
    end
%% PGE minimization approach 
else   
    ite   = 1;
    for lp=1:4
        for i=size(C,2):size(C_,1)
            B1      = A1*C_(i,:)';
            dG      = G_(i) - sum(B1.*g_A');

            ph2swp  = -1;
            if dG < -1e-8
                min_F  =  1e6;
                for j=1:size(C,2)
                    F = n_vec(j)/B1(j);
                    if F < min_F && F > 0.0
                        min_F  = F;
                        ph2swp = j;
                    end
                end
            end

            if ph2swp ~= -1
                swp     = 1;
                n_swap  = n_swap + 1;
                idA(ph2swp) = i;
                g_A(ph2swp)     = G_(i);

                for j=1:size(C,2)
                    A(:,ph2swp) = C_(i,:)';
                end

                A0      = A;
                A1      = inv(A0);   
                n_vec   = A\b';
            end

        end  
    end
    Gamma = (A'\g_A')';
    
    %Partitioning Gibbs Energy stage (PGE loops)
    for pge_ite=1:max_ite+1  
        xp = [];
        rot   = rotate_ghpp(C,Gamma);                                       % apply change of base to the pure endmembers
        for j=1:size(idA,2)
            if idA(j) > n_pp                                                % if the pseudocompound does not belong to a fake phase or a pure endmember
                sol             = fmincon(@(x) calc_G_fct(rot,x,gb,v,W,T,R,facpl),[x_(idA(j),1) x_(idA(j),2)],[],[],[],[],[1e-8 1e-8],[1-1e-8 1-1e-8],[],options);
                [G,phC]         = calc_GC_fct(zeros(1,size(C,2)),[sol(1) sol(2)],C,gb,v,W,T,R,facpl);
                G_(idA(j))      = G;
                C_(idA(j),:)    = phC(1,:);
                xx_(idA(j),:)   = [sol(1) sol(2)];
                color           = [color, min(0.0+(0.9/max_ite*3)*pge_ite,1)];
                x_              = [x_; [sol(1) sol(2)]];
                
                p               = x2p_fct([sol(1) sol(2)]);       
                if pge_ite > 1
                    xp          = [xp, [sol(1) sol(2)]];
                end
                [G,phC,xi]      = calc_GC_fct(rot,[sol(1) sol(2)],C,gb,v,W,T,R,facpl);
                xi_(idA(j),:)   = xi(1,:);
                p_(idA(j),:)    = p(1,:);
            end

        end
        min_path = [min_path; xp];

        % Not a generalized way to remove phases that converged to the same point
        for j=2:size(idA,2)
            if abs(G_(idA(j)) - G_(idA(j-1))) < 1e-5
                idA(j)      = [];
                n_vec(j-1)  = n_vec(j-1) + n_vec(j);
                n_vec(j)    = [];
                break
            end
        end

        % get number of active pure phases and solution phases
        PP = 0;
        SS = 0;
        for j=1:size(idA,2)
            if idA(j) <= n_pp   %if the phase is a pure phases
                PP = PP + 1;
            end
            if idA(j) > n_pp   %if the phase is a pure phases
                SS = SS + 1;
            end  
        end

        % build PGE newton-raphson system
        for i=1:n_C
            RHStrue(i) = b(i);
            for j=1:size(idA,2)
                if idA(j) <= n_pp   %if the phase is a pure phases
                    RHStrue(i) = RHStrue(i) - n_vec(j)*C_(idA(j),i); 
                end
                if idA(j) > n_pp
                    for k=1:n_em
                        RHStrue(i) = RHStrue(i) - n_vec(j)*p_(idA(j),k)*C(k,i)*facpl;
                    end
                end
            end
        end
        massR(:,pge_ite) = norm(RHStrue(:));
        for i=1:n_C
            RHS(i) = b(i);
            for j=1:size(idA,2)
                if idA(j) <= n_pp   %if the phase is a pure phases
                    RHS(i) = RHS(i) - n_vec(j)*C_(idA(j),i); 
                end
                if idA(j) > n_pp
                    for k=1:n_em
                        RHS(i) = RHS(i) - n_vec(j)*xi_(idA(j),k)*C(k,i)*facpl;
                    end
                end
            end
        end
        for k=1:SS
            RHS(n_C+k)      = 1 - sum( xi_(idA(k+PP),:) );
        end
        for k=1:PP
            RHS(n_C+SS+k)   = G_(idA(k)) - C_(idA(k),:)*Gamma';
        end

        LHS = zeros(n_C+PP+SS);
        % Fill the xi formulation part of the matrix
        for k=1:n_C
            for j=1:n_C
                for i=1:size(idA,2)
                    if idA(i) > n_pp
                        for x=1:n_em
                            LHS(k,j) = LHS(k,j) + C(x,j)*facpl*C(x,k)*facpl*xi_(idA(i),x)*n_vec(i);
                        end
                    end
                end
            end
        end

        % sum xi = 1
        for k=1:SS
            for j=1:n_C
                for i=1:n_em
                    LHS(n_C+k,j) = LHS(n_C+k,j) + C(i,j)*facpl*xi_(idA(k+PP),i);
                    LHS(j,n_C+k) = LHS(n_C+k,j) + C(i,j)*facpl*xi_(idA(k+PP),i);
                end
            end
        end

        % pure phase part of the matrix
        for k=1:PP
            for j=1:n_C
                LHS(n_C+SS+k,j) = C_(idA(k),j);
                LHS(j,n_C+SS+k) = C_(idA(k),j);
            end
        end

        %solve the system
        PGEsol = LHS\RHS';

        %update solution
        for j=1:n_C
            Gamma(j) = Gamma(j) + PGEsol(j);
        end
        for j=1:SS
            n_vec(j+PP) = n_vec(j+PP) + PGEsol(j+n_C);
        end
        for j=1:PP
            n_vec(j) = n_vec(j) + PGEsol(j+n_C+SS);
        end
        
        % calculate residual with respect to the reference Gamma
        Residual(pge_ite) = norm(Gamma_sol-Gamma);
    end
end


%% Plot Gibbs surface
rot = rotate_ghpp(C,Gamma);
ii  = 1; 
it  = 1;
for i=0:0.0025:1
    jj = 1;
    for j=0:0.0025:1
        x           = [i j];
        sf          = sf_fct(x);
        sf_ok       = check_sf(sf);
        if sf_ok == 1
            G(ii,jj)    = calc_G_fct(rot,x,gb,v,W,T,R,facpl);  
            if (~isnan(G(ii,jj)))
                GT(1,it)    = 5+G(ii,jj);  
                XT(1,it)    = i;
                YT(1,it)    = j;
                it          = it + 1;
            end
        else
            G(ii,jj) = NaN;
        end
        X(ii,jj)    = i;
        Y(ii,jj)    = j;
        jj          = jj + 1;
    end
    ii = ii + 1;
end

figure(1)
subplot(121)
[cA,cB] = c2t(X,Y);
cA(find(isnan(G))) = NaN;
cB(find(isnan(G))) = NaN;
h       = pcolor(cA,cB,G);
hold on
v       = [0.0:0.25:2.5];
vv      = [-0.5:0.25:0.0];
vvv     = [0 0];
contour(cA,cB,G,v,'w--')
contour(cA,cB,G,vv,'w--')
contour(cA,cB,G,vvv,'w-')
set(h, 'EdgeColor', 'none');

for i=size(C,2):size(x_,1)
    [pA,pB] = c2t(x_(i,1),x_(i,2));
    plot(pA,pB,'.','Color', [color(i) color(i) color(i)],'MarkerSize',12)
end

c               = colorbar;
c.Label.String  = '\DeltaG';
daspect([1 1 1])
xlabel('X','fontweight','bold') 
ylabel('Y','fontweight','bold')

str = '+qtz +sill';
text(0.75, 0.8, str,'FontSize',11)
str = 'An';
text(0.5, 0.9, str,  'fontweight','bold','FontSize',11)
str = 'Or';
text(1.0, -0.05, str,    'fontweight','bold','FontSize',11)
str = 'Ab';
text(-0.1, -0.05, str,  'fontweight','bold','FontSize',11)
set(gca,'Visible','off')
hold off


figure(1)
% Gamma residual
Gamma_lp  = [3.61075225756705,1.99399830509913,0.154754011502646,0.796900504123559,0.339236902549688,0.114164049993410,0.00230974725952924,0.0539473167380940,0.0256808398330383,0.0115332774624007,0.00482469098673758,0.00148631941233943,0.000186120276705384,0.000645148625750794,0.000226467950795291,1.84883537479303e-05,1.84883537479303e-05,1.84883537479303e-05,1.84883537479303e-05,1.84883537479303e-05,1.84883537479303e-05,1.84883537479303e-05,1.84883537479303e-05,1.84883537479303e-05];
Gamma_PGE = [3.61075737154559,0.901386432406202,0.422509908270587,0.194632500975321,0.0889871537077892,0.0405498980931909,0.0184517459340025,0.00838997109461454,0.00381454625185100,0.00173354353453510,0.000788094857296457,0.000358239174337608,0.000163091009525831,7.52573999074663e-05,3.37601238260706e-05,1.53689699434672e-05,1.47145102125727e-05,9.24311698934173e-06,3.72353036029478e-06,3.74236734960299e-06,2.95625618851991e-07,4.44776890384541e-06,4.51650338843106e-06,3.66685703268728e-06,1.62118493169185e-07];
subplot(122)
plot(log10(Gamma_lp),'ks--','MarkerFaceColor','w')
hold on
plot(log10(Gamma_PGE),'ko-','MarkerFaceColor','w')
title('Convergence profile: Linear Programming vs PGE')
xlabel('Iteration number','fontweight','bold') 
ylabel('log10(norm of \Gamma residual)','fontweight','bold')
legend('LP','PGE')
hold off

pl1LP = [0.0102948150000000,0.892833570000000;0.0106091450000000,0.818086984000000;0.0109390720000000,0.845518634000000;0.0104948310000000,0.836404278000000;0.0106124640000000,0.842405981000000;0.0106757750000000,0.845224224000000;0.0107084790000000,0.846594949000000;0.0107255430000000,0.847273172000000;0.0107164810000000,0.846934887000000;0.0107119950000000,0.846765581000000;0.0107097730000000,0.846685745000000;0.0107087710000000,0.846645226000000;0.0107082210000000,0.846625269000000;0.0107085000000000,0.846635490000000;0.0107082840000000,0.846628437000000;0.0107082840000000,0.846628437000000;0.0107082840000000,0.846628437000000;0.0107082840000000,0.846628437000000;0.0107082840000000,0.846628437000000;0.0107082840000000,0.846628437000000;0.0107082840000000,0.846628437000000;0.0107082840000000,0.846628437000000;0.0107082840000000,0.846628437000000];
pl2LP = [0.597721605000000,0.00406070100000000;0.338867614000000,0.0157856370000000;0.433204155000000,0.00965875600000000;0.385447881000000,0.0126569690000000;0.409210329000000,0.0111478510000000;0.421168562000000,0.0104460230000000;0.427174585000000,0.0101075070000000;0.430203298000000,0.00994023300000000;0.428667012000000,0.0100250430000000;0.427900333000000,0.0100675980000000;0.427536643000000,0.0100878350000000;0.427355703000000,0.0100979530000000;0.427309873000000,0.0101004860000000;0.427287266000000,0.0101016740000000;0.427276722000000,0.0101023560000000;0.427276433000000,0.0101023310000000;0.427270401000000,0.0101027200000000;0.427270401000000,0.0101027200000000;0.427270401000000,0.0101027200000000;0.427270401000000,0.0101027200000000;0.427270401000000,0.0101027200000000;0.427270401000000,0.0101027200000000;0.427270401000000,0.0101027200000000];

pl1PGE = [0.0102948130000000,0.892833603000000;0.0106672470000000,0.848976171000000;0.0106854130000000,0.847623324000000;0.0106977130000000,0.847091884000000;0.0107035930000000,0.846848747000000;0.0107062610000000,0.846733772000000;0.0107074820000000,0.846677480000000;0.0107080350000000,0.846651222000000;0.0107081600000000,0.846638562000000;0.0107082080000000,0.846633986000000;0.0107082950000000,0.846630286000000;0.0107082000000000,0.846629252000000;0.0107084050000000,0.846627605000000;0.0107082600000000,0.846628637000000;0.0107082720000000,0.846628270000000;0.0107083660000000,0.846627052000000;0.0107082600000000,0.846628307000000;0.0107082480000000,0.846628407000000;0.0107083640000000,0.846626844000000;0.0107082550000000,0.846628338000000;0.0107082020000000,0.846628224000000;0.0107082440000000,0.846628357000000;0.0107082610000000,0.846628308000000;0.0107082900000000,0.846627145000000;0.0107082540000000,0.846628194000000];
pl2PGE = [0.597720254000000,0.00406074400000000;0.433581737000000,0.00980383900000000;0.429703249000000,0.00999095900000000;0.428402859000000,0.0100504960000000;0.427827578000000,0.0100766910000000;0.427550937000000,0.0100893780000000;0.427411647000000,0.0100958690000000;0.427344088000000,0.0100990210000000;0.427308475000000,0.0101008070000000;0.427291645000000,0.0101016170000000;0.427283102000000,0.0101020470000000;0.427276957000000,0.0101023730000000;0.427276050000000,0.0101024110000000;0.427276982000000,0.0101023110000000;0.427275225000000,0.0101024300000000;0.427273305000000,0.0101025400000000;0.427274440000000,0.0101024920000000;0.427274885000000,0.0101024660000000;0.427275020000000,0.0101024720000000;0.427274976000000,0.0101024250000000;0.427274013000000,0.0101025140000000;0.427274248000000,0.0101024830000000;0.427274369000000,0.0101024600000000;0.427275532000000,0.0101024410000000;0.427274724000000,0.0101024490000000];

figure(2)
subplot(121)
plot(pl2LP(:,1),'k--s','MarkerFaceColor','w')
hold on
plot(pl2PGE(:,1),'k-o','MarkerFaceColor','w')
subplot(122)
plot(pl1LP(:,2),'k--s','MarkerFaceColor','w')
hold on
plot(pl1PGE(:,2),'k-o','MarkerFaceColor','w')

%% Function definitions

%Generate the pseudocompounds
function [G_,C_,x_] = gen_PC_fct(Gamma,C,gb,v,W,T,R,facpl,gbsil,Csil,facsil,gbqtz,Cqtz,facqtz)
    G_ = [];
    C_ = [];
    x_ = [];
    
    %add starting penalty oxides
    for i=1:size(C,2)
        G_ = [G_,1e6];
        Cfake = zeros(1,5);
        Cfake(1,i) = 1;
        C_ = [C_;Cfake];
        x_ = [x_;[-1 -1]];
    end  
    
    %add pure endmembers
    for i=1:size(gb,2)
        G_ = [G_,gb(i)*facpl];
        C_ = [C_;C(i,:)*facpl];
    end
    x_ = [x_;[0 0]];
    x_ = [x_;[1 0]];
    x_ = [x_;[0 1]];
        
    G_ = [G_,gbqtz*facqtz];
    C_ = [C_;Cqtz*facqtz];
    x_ = [x_;[-1 -1]];    
    G_ = [G_,gbsil*facsil];
    C_ = [C_;Csil*facsil];
    x_ = [x_;[-1 -1]];
    
    %Assemble pseudocompounds
    min_v   = 1e-6;                                                         %slight shift from the bounds
    step    = 0.249;                                                        %discretization step
    cor     = 1/step*min_v;

    i = 1;
    for x = min_v:step-cor:1.0
        for y = min_v:step-cor:1.0
            sf      = sf_fct([x y]);
            sf_ok   = check_sf(sf);
            if sf_ok == 1
                PC(i,:)    = [x y];
                i = i+1;   
            end

        end
    end

    for j =1:size(PC,1)
        [G,phC] = calc_GC_fct(Gamma,[PC(j,1) PC(j,2)],C,gb,v,W,T,R,facpl);
        G_      = [G_,G];
        C_      = [C_;phC];
        x_      = [x_;[PC(j,1) PC(j,2)]];
    end

end

%hyper-plane rotation, change of base
function rot = rotate_ghpp(C,Gamma)
    rot(1) = sum(C(1,:).*Gamma);
    rot(2) = sum(C(2,:).*Gamma);
    rot(3) = sum(C(3,:).*Gamma);
end

% calculation Gibbs energy of pl4T
function G      = calc_G_fct(rot,x,gb,v,W,T,R,facpl)
    p           = x2p_fct(x);
    Gex         = Gex_fct(p,v,W);
    sf          = sf_fct(x);
    mu          = mu_fct(rot,T,R,sf,gb,Gex);
    G           = Gpl_fct(mu,p,facpl);
end

% calculation Gibbs energy of pl4T + getting back composition and xi
function [G,phC,xi]  = calc_GC_fct(rot,x,C,gb,v,W,T,R,facpl)
    p           = x2p_fct(x);
    Gex         = Gex_fct(p,v,W);
    sf          = sf_fct(x);
    mu          = mu_fct(rot,T,R,sf,gb,Gex);
	xi          = xi_fct(mu,p,T,R);
    G           = Gpl_fct(mu,p,facpl);
    phC         = get_phC_fct(p,C,facpl);
end

% get phase composition
function phC = get_phC_fct(p,C,facpl)
    phC = p*C*facpl;
end

% composition variables to endmember fractions
function p = x2p_fct(x)
    p(1) = - x(1) - x(2) + 1.0;
    p(2) =  x(1);
    p(3) =  x(2);
end

% G excess calculation
function Gex = Gex_fct(p,v,W)
    E        = eye(size(p,2));
    sum_v    = sum(p.*v);
    mat_phi  = p.*v/sum_v;
 
    for i=1:size(p,2)
        Gex(i)  = 0.0;
        it      = 1;
        for j=1:size(p,2)-1
            for k=j+1:size(p,2)
                Gex(i) = Gex(i) - (E(i,j) - mat_phi(j))*(E(i,k) - mat_phi(k))*(W(it)*2.0*v(i)/(v(j)+v(k)));
                it = it + 1;
            end
        end
    end
    
end

% Site-fractions calculation
function sf = sf_fct(x)
    sf(1)           = -x(1) - x(2) + 1.0;
    sf(2)           =  x(1);
    sf(3)           =  x(2);
    sf(4)           =  0.25*x(1) + 0.25;
    sf(5)           =  0.75 - 0.25*x(1);
end

% Site-fractions check
function sf_ok = check_sf(sf)
    sf_ok = 1;
    if any(sf < 0.0)
        sf_ok = 0;
    end
end

% Chemical potential calculation
function mu = mu_fct(rot,T,R,sf,gb,Gex)
	mu(1)          = R*T*real(log(1.7548*sf(1)*power(sf(4), 0.25)*power(sf(5), 0.75))) 	+ gb(1) - rot(1) + Gex(1);
	mu(2)          = R*T*real(log(2.0*sf(2)*sqrt(sf(4))*sqrt(sf(5))))                   + gb(2) - rot(2) + Gex(2);
	mu(3)          = R*T*real(log(1.7548*sf(3)*power(sf(4), 0.25)*power(sf(5), 0.75))) 	+ gb(3) - rot(3) + Gex(3);
end

%Gibbs energy calculation
function Gpl = Gpl_fct(mu,p,facpl)
    Gpl = sum(p.*mu)*facpl;
end

%xi calculation
function xip = sum_pen(mu,p,T,R)
    xip = sum(exp(-mu./(R*T)).*p);
end

%xi calculation
function xi = xi_fct(mu,p,T,R)
    xi = exp(-mu./(R*T)).*p;
end

%cart2term
function [a,b] = c2t(X,Y)
    a      = 0.5*(2*X+Y);
    b      = (sqrt(3)/2)*Y;
end
