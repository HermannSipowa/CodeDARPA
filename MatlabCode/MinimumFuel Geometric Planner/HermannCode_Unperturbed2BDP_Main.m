%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hermann Kaptui Sipowa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
clear
close all
start_up
format long e
clc
global mu_Earth tspan M T_normalized mu_dot ...
    Rrho Tc Num_Agent a_chief e_chief DynamicsModel ...
    Period % iteration InitalCondition

c1  = rgb('DarkGreen');
c2  = rgb('Gold');
c3  = rgb('Lime');
c4  = rgb('DarkOrange');
c5  = rgb('DarkBlue');
c6  = rgb('Red');
c7  = rgb('Purple');
c8  = rgb('Bisque');
c9  = rgb('Orange');
c10 = rgb('DarkGray');
c11 = rgb('Teal');
mm = 6;
addpath('../casadi-osx-matlabR2015a-v3.5.5')
import casadi.*
CasADiopts = struct('main', true, 'mex', true);

%%
%------------------- Specify the dynamics Model ------------------------  
% Chose from: 
% - CWode = Clohessy–Wiltshire, 
% - THode = Tschauner-Hempel , 
% - NLode = Full Nonlinear Model 
DynamicsModel = 'CWode';
SimTime = 0.75;

%% Calculating the boundaries conditions
%=======================================%
mu_Earth   = 3.986004415E5; % Earth gravitational parameter
Num_Agent  = 2; % Number of agents in the system
Num_Agent1 = 1; % Number of agents in the first final-orbit
Num_Agent2 = 1; % Number of agents in the second final-orbit

% Specifying the chief's orbit
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
a_chief      = 1.42e4;   % Semi-major axis in Km
e_chief      = 0.0;      % Eccentricity
inc_chief    = 50;       % Inclination in deg0
BigOmg_chief = 10;       % RAAN in deg
LitOmg_chief = 10;       % AOP in deg
M_chief1     = 0;        % Initial mean anomaly

%% Integrating the chief's trajectory
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
Period  = (2*pi)*sqrt(a_chief^3/mu_Earth);
IntTime = SimTime*Period;
tspan   = linspace(0,IntTime,1e4);
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-17);
mu_dot  = sqrt(mu_Earth/a_chief^3); % Chief's mean motion
Rc  = a_chief*eye(3); Rrho = eye(3); Tc = Period; T_normalized = tspan/Tc;
inc_chief = deg2rad(inc_chief); BigOmg_chief = deg2rad(BigOmg_chief);
LitOmg_chief = deg2rad(LitOmg_chief); M_chief1 = deg2rad(M_chief1);
COE1 = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief1];
[Position_target1,Velocity_target1] = COEstoRV(COE1,mu_Earth);
X0_Chief1  = [Position_target1; Velocity_target1];
[~, Xnom] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),tspan,X0_Chief1,options);

Ttraj = linspace(0,Period,1e4); X0_Chief2  = Xnom(end,:)'; 
[~, Xnom0] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),Ttraj,X0_Chief1,options);
[~, Xnomf] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),Ttraj,X0_Chief2,options);



%% ------- Create MX files to interpolate the chief's states (Xnom) -------
for i = 1
    if exist('Chief_x','file') == 3
        delete Chief_x.c Chief_x.mexmaci64
    end
    if exist('Chief_y','file') == 3
        delete Chief_y.c Chief_y.mexmaci64
    end
    if exist('Chief_z','file') == 3
        delete Chief_z.c Chief_z.mexmaci64
    end
    if exist('Chief_dx','file') == 3
        delete Chief_dx.c Chief_dx.mexmaci64
    end
    if exist('Chief_dy','file') == 3
        delete Chief_dy.c Chief_dy.mexmaci64
    end
    if exist('Chief_dz','file') == 3
        delete Chief_dz.c Chief_dz.mexmaci64
    end
    
    Chief_x  = casadi.interpolant('Chief_x','bspline',{T_normalized},Xnom(:,1)');
    Chief_y  = casadi.interpolant('Chief_y','bspline',{T_normalized},Xnom(:,2)');
    Chief_z  = casadi.interpolant('Chief_z','bspline',{T_normalized},Xnom(:,3)');
    Chief_dx = casadi.interpolant('Chief_dx','bspline',{T_normalized},Xnom(:,4)');
    Chief_dy = casadi.interpolant('Chief_dy','bspline',{T_normalized},Xnom(:,5)');
    Chief_dz = casadi.interpolant('Chief_dz','bspline',{T_normalized},Xnom(:,6)');
    
    Chief_x.generate('Chief_x.c',CasADiopts);
    Chief_y.generate('Chief_y.c',CasADiopts);
    Chief_z.generate('Chief_z.c',CasADiopts);
    Chief_dx.generate('Chief_dx.c',CasADiopts);
    Chief_dy.generate('Chief_dy.c',CasADiopts);
    Chief_dz.generate('Chief_dz.c',CasADiopts);
    mex Chief_x.c -largeArrayDims
    mex Chief_y.c -largeArrayDims
    mex Chief_z.c -largeArrayDims
    mex Chief_dx.c -largeArrayDims
    mex Chief_dy.c -largeArrayDims
    mex Chief_dz.c -largeArrayDims
    clear Chief_x Chief_y Chief_z Chief_dx Chief_dy Chief_dz
    
    
    if strcmp(DynamicsModel,'NLode')
        if exist('Chief_x0','file') == 3
            delete Chief_x0.c Chief_x0.mexmaci64
        end
        if exist('Chief_y0','file') == 3
            delete Chief_y0.c Chief_y0.mexmaci64
        end
        if exist('Chief_z0','file') == 3
            delete Chief_z0.c Chief_z0.mexmaci64
        end
        if exist('Chief_dx0','file') == 3
            delete Chief_dx0.c Chief_dx0.mexmaci64
        end
        if exist('Chief_dy0','file') == 3
            delete Chief_dy0.c Chief_dy0.mexmaci64
        end
        if exist('Chief_dz0','file') == 3
            delete Chief_dz0.c Chief_dz0.mexmaci64
        end
        TT = Ttraj/Tc;
        Chief_x0  = casadi.interpolant('Chief_x0','bspline',{TT},Xnom0(:,1)');
        Chief_y0  = casadi.interpolant('Chief_y0','bspline',{TT},Xnom0(:,2)');
        Chief_z0  = casadi.interpolant('Chief_z0','bspline',{TT},Xnom0(:,3)');
        Chief_dx0 = casadi.interpolant('Chief_dx0','bspline',{TT},Xnom0(:,4)');
        Chief_dy0 = casadi.interpolant('Chief_dy0','bspline',{TT},Xnom0(:,5)');
        Chief_dz0 = casadi.interpolant('Chief_dz0','bspline',{TT},Xnom0(:,6)');
    
        Chief_x0.generate('Chief_x0.c',CasADiopts);
        Chief_y0.generate('Chief_y0.c',CasADiopts);
        Chief_z0.generate('Chief_z0.c',CasADiopts);
        Chief_dx0.generate('Chief_dx0.c',CasADiopts);
        Chief_dy0.generate('Chief_dy0.c',CasADiopts);
        Chief_dz0.generate('Chief_dz0.c',CasADiopts);
        mex Chief_x0.c -largeArrayDims
        mex Chief_y0.c -largeArrayDims
        mex Chief_z0.c -largeArrayDims
        mex Chief_dx0.c -largeArrayDims
        mex Chief_dy0.c -largeArrayDims
        mex Chief_dz0.c -largeArrayDims
        clear Chief_x0 Chief_y0 Chief_z0 Chief_dx0 Chief_dy0 Chief_dz0
    
    
        if exist('Chief_xf','file') == 3
            delete Chief_xf.c Chief_xf.mexmaci64
        end
        if exist('Chief_yf','file') == 3
            delete Chief_yf.c Chief_yf.mexmaci64
        end
        if exist('Chief_zf','file') == 3
            delete Chief_zf.c Chief_zf.mexmaci64
        end
        if exist('Chief_dxf','file') == 3
            delete Chief_dxf.c Chief_dxf.mexmaci64
        end
        if exist('Chief_dyf','file') == 3
            delete Chief_dyf.c Chief_dyf.mexmaci64
        end
        if exist('Chief_dzf','file') == 3
            delete Chief_dzf.c Chief_dzf.mexmaci64
        end
    
        Chief_xf  = casadi.interpolant('Chief_xf','bspline',{TT},Xnomf(:,1)');
        Chief_yf  = casadi.interpolant('Chief_yf','bspline',{TT},Xnomf(:,2)');
        Chief_zf  = casadi.interpolant('Chief_zf','bspline',{TT},Xnomf(:,3)');
        Chief_dxf = casadi.interpolant('Chief_dxf','bspline',{TT},Xnomf(:,4)');
        Chief_dyf = casadi.interpolant('Chief_dyf','bspline',{TT},Xnomf(:,5)');
        Chief_dzf = casadi.interpolant('Chief_dzf','bspline',{TT},Xnomf(:,6)');
    
        Chief_xf.generate('Chief_xf.c',CasADiopts);
        Chief_yf.generate('Chief_yf.c',CasADiopts);
        Chief_zf.generate('Chief_zf.c',CasADiopts);
        Chief_dxf.generate('Chief_dxf.c',CasADiopts);
        Chief_dyf.generate('Chief_dyf.c',CasADiopts);
        Chief_dzf.generate('Chief_dzf.c',CasADiopts);
        mex Chief_xf.c -largeArrayDims
        mex Chief_yf.c -largeArrayDims
        mex Chief_zf.c -largeArrayDims
        mex Chief_dxf.c -largeArrayDims
        mex Chief_dyf.c -largeArrayDims
        mex Chief_dzf.c -largeArrayDims
        clear Chief_xf Chief_yf Chief_zf Chief_dxf Chief_dyf Chief_dzf
    end
end 
% -------------------------------------------------------------------------

%% --------------- Setting the deputy initial conditions ------------------
for i = 1
    r0 = Xnom0(1,1:3)'; v0 = Xnom0(1,4:6)';
    COE0 = RVtoCOEs(r0,v0,mu_Earth); f1 = COE0(6);
    q1_chief = e_chief*cos(LitOmg_chief); q2_chief = e_chief*sin(LitOmg_chief);
    theta_chief1 = f1+LitOmg_chief;
    OE_chief1 = [a_chief, theta_chief1, inc_chief, q1_chief, q2_chief, BigOmg_chief];
    AMap1 = ForwardMapping(OE_chief1, mu_Earth); % Linear mapping matrix
    
    rf = Xnomf(1,1:3)'; vf = Xnomf(1,4:6)';
    COEf = RVtoCOEs(rf,vf,mu_Earth); f2 = COEf(6);
    theta_chief2 = f2+LitOmg_chief;
    OE_chief2 = [a_chief, theta_chief2, inc_chief, q1_chief, q2_chief, BigOmg_chief];
    AMap2 = ForwardMapping(OE_chief2, mu_Earth); % Linear mapping matrix
    
    % Specify the final relative-orbital elements of the deputies in both orbits
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    dela = 0;
    dele1 = 1/(6*a_chief);
    deli1 = -1/(3*a_chief);
    delLitOmg = -2*pi*1E-5;
    delBigOmg = 0;
    delM = pi*1E-6;
    delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_1 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
    
    
    dele1 = -1/(8*a_chief);
    deli1 = -1/(20*a_chief);
    delLitOmg = -5*pi*1E-7;
    delBigOmg = pi*1E-6;
    delM = 0;
    delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_2 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
    
    
    dele1 = 1/(6*a_chief);
    deli1 = 1/(2*a_chief);
    delLitOmg = 0;
    delBigOmg = 0;
    delM = 0;
    delq1_1 = dele1*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele1*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_3 = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
    
    dele2 = 1/(8*a_chief);
    deli2 = 2/(10*a_chief);
    delBigOmg = -pi*1E-5;
    delM = pi*1E-6;
    delLitOmg = pi*1E-5;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    
    delq1_2 = dele2*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_2 = dele2*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    delCOE_4 = [dela, deltheta, deli2, delq1_2, delq2_2, delBigOmg]';
    
    % Integrating the deputy initial and final trajectories
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    X01 = AMap1*delCOE_1; Xi1 = X01;
    X02 = AMap1*delCOE_2; Xi2 = X02;
    Xf1 = AMap2*delCOE_3;
    Xf2 = AMap2*delCOE_4;
    
end 
% -------------------------------------------------------------------------

%% ------------------------ AutoDiff using CasADi -------------------------
for j = 1
    XChief= SX.sym('XC',6); t = SX.sym('t'); Penalty = SX.sym('Penalty');
    rho1  = SX.sym('rho1',6); rho2    = SX.sym('rho2',6);
    drho1 = SX.sym('rho1',6); drho2   = SX.sym('rho2',6);
    Xaug  = [rho1;rho2];      dXaug   = [drho1;drho2];
    M = 3; N = 6; drift = SX(N*Num_Agent,1);
    
    if exist('L','file') == 3
        delete L.c L.mexmaci64
    end
    if exist('dLdx','file') == 3
        delete dLdx.c dLdx.mexmaci64
    end
    if exist('fdrift','file') == 3
        delete fdrift.c fdrift.mexmaci64
    end
    
    for i = 1:Num_Agent
        if i == 1
            % [Fc F], note: the F_bar matrix in paper
            F = eye(N); F_Aug = F;
            % penalty matrix (D matrix in the paper)
            D = diag([Penalty*ones(1,N-M) ones(1,M)]); D_Aug = D;
        else
            F_Aug = blkdiag(F_Aug,F);
            D_Aug = blkdiag(D_Aug,D);
        end
    end
    
    % Defining the drift vector field
    for i = 1:Num_Agent
        idx = 1+N*(i-1):N*i;
        if strcmp(DynamicsModel,'THode')
            drift(idx,1) = NonDimentionalized_THodeCasADi(t,Xaug(idx),XChief);
        elseif strcmp(DynamicsModel,'CWode')
            drift(idx,1) = NonDimentionalized_CWodeCasADi(t,Xaug(idx),XChief);
        elseif strcmp(DynamicsModel,'NLode')
            drift(idx,1) = NonDimentionalized_NLodeCasADi(t,Xaug(idx),XChief);
        end
    end
    
    % Defining the metric, before adding the state constraint barrier functions
    H = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);
    
    % ------------------- state constraint barrier function ------------------
    B = [];
    % b = 0; % no state contraintst
    % B = [B b]; % attach the barrier function to B
    
    % -------------------- Metric and curve Length ----------------------------
    % the metric with state constraint barier functions
    G = (sum(B)+1)*H;
    % actuated curve length
    L = Function('L',{Xaug,dXaug,Penalty,t,XChief},...
        {(dXaug - drift).' * G * (dXaug - drift)},...
        {'X','dX','k','t','XChief'},...
        {'CurveLength'});
    dLdx = L.jacobian();
    fdrift = Function('fdrift', {t,Xaug,XChief},...
        {drift});
    
    % Rewrite the functions in c files
    L.generate('L.c',CasADiopts);
    dLdx.generate('dLdx.c',CasADiopts);
    fdrift.generate('fdrift.c',CasADiopts);
    
    % Generating Mex files from the c files
    mex L.c -largeArrayDims
    mex dLdx.c -largeArrayDims
    mex fdrift.c -largeArrayDims
    
    clear L dLdx fdrift
end
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% ##################### Geometric Motion Planning #########################
% -------------------------------------------------------------------------
%% ------- Nondimentionalized the inital condition of the problem ---------
for j = 1
    
    if strcmp(DynamicsModel,'THode')
        h = sqrt(mu_Earth*a_chief*(1-e_chief^2));
        K_f = 1+e_chief*cos(f1);
        conts1 = mu_Earth*K_f/(h^(3/2));
        conts2 = -mu_Earth*e_chief*sin(f1)/(h^(3/2));
        conts3 = 1/conts1;
        A_map1 = [conts1*eye(3)       zeros(3);
                  conts2*eye(3)  conts3*eye(3)];
        x_bar01 = A_map1*X01;
        x_bar02 = A_map1*X02;
        X0  = [x_bar01; x_bar02];
        tspan1 = 2*pi*T_normalized;
        h      = sqrt(mu_Earth*a_chief*(1-e_chief^2));
        K_f    = 1+e_chief*cos(f2);
        conts1 = mu_Earth*K_f/(h^(3/2));
        conts2 = -mu_Earth*e_chief*sin(f2)/(h^(3/2)); 
        conts3 = 1/conts1;
        A_map2 = [conts1*eye(3)       zeros(3);
                  conts2*eye(3)  conts3*eye(3)];
        x_bar11 = A_map2*Xf1; 
        x_bar21 = A_map2*Xf2; 
        Xf = [x_bar11; x_bar21];
        tspan2 = linspace(f2,(f2+2*pi),length(T_normalized));
        tic
        [~, X_rel_0] = ode113(@(t,x)NonDimentionalized_THode(t,x),tspan1,X0,options);
        [~, X_rel_1] = ode113(@(t,x)NonDimentionalized_THode(t,x),tspan2,Xf,options);
        toc
        X_THpos0 = nan(length(X_rel_0),6);
        X_THpos1 = nan(length(X_rel_0),6);
        X_THpos2 = nan(length(X_rel_0),6);
        X_THpos3 = nan(length(X_rel_0),6);
        for i = 1:length(X_rel_0)
            f = tspan1(i);
            h = sqrt(mu_Earth*a_chief*(1-e_chief^2));
            K_f = 1+e_chief*cos(f);
            conts1 = h^(3/2)/(mu_Earth*K_f);
            conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
            conts3 = 1/conts1;
            A_inv1 = [conts1*eye(3)       zeros(3);
                      conts2*eye(3)  conts3*eye(3)];
            X_THpos0(i,:) = (A_inv1*X_rel_0(i,1:6).').';
            X_THpos1(i,:) = (A_inv1*X_rel_0(i,7:12).').';
            
            f = tspan2(i);
            h = sqrt(mu_Earth*a_chief*(1-e_chief^2));
            K_f = 1+e_chief*cos(f);
            conts1 = h^(3/2)/(mu_Earth*K_f);
            conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
            conts3 = 1/conts1;
            A_inv2 = [conts1*eye(3)       zeros(3);
                      conts2*eye(3)  conts3*eye(3)];
            X_THpos2(i,:) = (A_inv2*X_rel_1(i,1:6).').';
            X_THpos3(i,:) = (A_inv2*X_rel_1(i,7:12).').';
        end
        X_THpos01 = [X_THpos0 X_THpos1];
        X_THpos11 = [X_THpos2 X_THpos3];
    
    elseif strcmp(DynamicsModel,'CWode')
         X0 = [X01; X02]/Tc;
         Xf = [Xf1; Xf2]/Tc;
         TT = Ttraj/Tc;
         [~, X_rel_0] = ode113(@(t,x)NonDimentionalized_CWode(t,x),TT,X0);
         [~, X_rel_1] = ode113(@(t,x)NonDimentionalized_CWode(t,x),TT,Xf);
         X_THpos01 = Tc*X_rel_0;
         X_THpos11 = Tc*X_rel_1;
             
    elseif strcmp(DynamicsModel,'NLode')
        X0 = [X01 X02];
        Xf = [Xf1 Xf2];
        tic
        options = odeset('RelTol',2.22045e-9,'AbsTol',2.22045e-30);
        X_THpos01 = []; X_THpos11 = []; X_rel_0 = []; X_rel_1 = [];
        for i = 1:2
            X0integral = [X01,X02];
            Xfintegral = [Xf1,Xf2];
            TT = Ttraj/Tc;
            [~, X_rel_01] = ode113(@(t,x)NonDimentionalized_NLode0(t,x),TT,X0(:,i),options);
            [~, X_rel_11] = ode113(@(t,x)NonDimentionalized_NLodef(t,x),TT,Xf(:,i),options);
            
            X_rel_0   = [X_rel_0, X_rel_01];  X_rel_1 = [X_rel_1, X_rel_11];
            X_THpos01 = [X_THpos01, (Rrho*X_rel_01(:,1:3).').', (Rrho*X_rel_01(:,4:6).').'];
            X_THpos11 = [X_THpos11, (Rrho*X_rel_11(:,1:3).').', (Rrho*X_rel_11(:,4:6).').'];
        end
        toc
        X0 = X0(:); Xf = Xf(:);
    end
    
end
% -------------------------------------------------------------------------

%% --------------------------- AGHF parameters ----------------------------
for ll = 1
    clc
    smax = 5e7; % Need tpo find an analytic way to solve for this value
    % # of grids in t
    tgrids = 250;
    % # of grids in s
    sgrids = 250;
    % # of grids of integration for final trajectory extraction
    intgrids = 1e5;
    % # of inputs
    M = 3;  N  = length(X0);
    % penalty value (the lambda in paper)
    Penalty = 5e5;
    % tolerance of the integrator
    % opts = odeset('AbsTol',1e-14);
    opts = odeset('RelTol',2.22045e-8,'AbsTol',2.22045e-20);
    
    % Setting the boundary condtions and the intial condition
    tmax = smax; xpoints = tgrids; xpoints_new = intgrids; tpoints = sgrids;
    m = 0; t = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of s interval in log scale
    
    if strcmp(DynamicsModel,'THode')
        T    = SimTime*wrapTo2Pi(f2)-f1; % motion duration
        x    = linspace(f1,SimTime*wrapTo2Pi(f2),xpoints); % discretization of the curve in time
        xnew = linspace(f1,SimTime*wrapTo2Pi(f2),intgrids);
        
    elseif strcmp(DynamicsModel,'CWode')
        T    = SimTime; % motion duration
        x    = linspace(0,T,xpoints); % discretization of the curve in time
        xnew = linspace(0,T,intgrids);
        
    elseif strcmp(DynamicsModel,'NLode')
        T    = SimTime; % motion duration
        x    = linspace(0,T,xpoints);
        xnew = linspace(0,T,intgrids);
        
    end
    
end
% -------------------------------------------------------------------------

%% ------------------------- Solving the AGHF -----------------------------
for ll = 1
    % solve for trajectory, see "AGHF" function file for implementation details
    close all
    tic;
    disp('solving PDE...');
    
    % Solve AGHF
    sol = pdepe(m,@(x,t,u,DuDx) mypdexpde(x,t,u,DuDx,Penalty,N),...
        @(x) mypdexic(x,X0,Xf,T),...
        @(xl,ul,xr,ur,t) mypdexbc(xl,ul,xr,ur,t,X0,Xf,N),...
        x,t);
    % The solution "sol" is of form sol(t,x,i)
    toc;
end
% -------------------------------------------------------------------------

%% ------------------ Calculating the action integral ---------------------
for ll = 1
    tic;
    disp('Computing the action integral...');
    % calculate the atuated curve length for each grid of s
    X_temp  = zeros(N,xpoints);
    dX_temp = zeros(N,xpoints);
    cost    = zeros(tpoints,1);
    
    n = 6;
    F = eye(n); F_Aug = blkdiag(F,F);
    % penalty matrix (D matrix in the paper)
    D = diag([Penalty*ones(1,n-M) ones(1,M)]); D_Aug = blkdiag(D,D);
    G = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);
    
    for j = 1:tpoints
        for kk = 1:N
            [X_temp(kk,:), dX_temp(kk,:)] = pdeval(m,x,sol(j,:,kk),x);
        end
        for i = 1:xpoints
            dX = dX_temp(:,i);
            for k = 1:Num_Agent
                idx = 1+n*(k-1):n*k;
                X = X_temp(idx,i);
                if strcmp(DynamicsModel,'THode')
                    f(idx,1) =  NonDimentionalized_THode(x(i),X);
                elseif strcmp(DynamicsModel,'CWode')
                    f(idx,1) = NonDimentionalized_CWode(x(i),X);
                elseif strcmp(DynamicsModel,'NLode')
                    f(idx,1) =  NonDimentionalized_NLode(x(i),X);
                end
            end
            cost(j) = cost(j) + (dX-f)'*G*(dX-f)*(T/(xpoints-1));
        end
    end
    toc;
end
% -------------------------------------------------------------------------

%% ----------------------- control extraction -----------------------------
for ll = 1
    disp('Extracting the required control...');
    % initialize controls and states
    Crtl = struct('u', cell(1, Num_Agent)); % actual control
    % use spline to interpolate state from steady state pde solution, p
    Pint   = zeros(N,xpoints_new);
    p      = zeros(N,xpoints_new);
    dpdx   = zeros(N,xpoints_new);
    dx_new = diff([xnew(3),xnew,xnew(xpoints_new-2)]);
    for i = 1:N
        Pint(i,:) = spline(x,sol(1,:,i),xnew);
        p(i,:)    = spline(x,sol(end,:,i),xnew);
        % Numerically computing the time derivative of the state
        dp_new    = diff([p(i,3),p(i,:),p(i,xpoints_new-2)]);
        dpdx(i,:) = (dp_new(1:end-1)./dx_new(1:end-1).*dx_new(2:end) ...
            + dp_new(2:end)./dx_new(2:end).*dx_new(1:end-1)) ...
            ./ (dx_new(2:end)+dx_new(1:end-1));
    end
    for i = 1 : length(xnew)
        if strcmp(DynamicsModel,'THode')
            f      = xnew(i);
            h      = sqrt(mu_Earth*a_chief*(1-e_chief^2));
            K_f    = 1+e_chief*cos(f);
            conts1 = h^(3/2)/(mu_Earth*K_f);
            conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
            conts3 = 1/conts1;
            A_inv  = [conts1*eye(3)       zeros(3);
                conts2*eye(3)  conts3*eye(3)];
    
            drift = NonDimentionalized_THode(xnew(i),p(:,i));
        elseif strcmp(DynamicsModel,'CWode')
            drift = NonDimentionalized_CWode(xnew(i),p(:,i));
        elseif strcmp(DynamicsModel,'NLode')
            d1 = NonDimentionalized_NLode(xnew(i),p(1:6,i));
            d2 = NonDimentionalized_NLode(xnew(i),p(7:12,i));
            drift = [d1;d2];
        end
        % get [Fc F], note: the F_bar matrix in paper
        Ffull = F;
        B = Ffull(:,M+1:end);
        % Extracting/approximate the required control
        for jj = 1:Num_Agent
            idx = 1+mm*(jj-1):mm*jj;
            Crtl(jj).u(:,i) = (B.'*B)\(B.') * ( dpdx(idx,i) - drift(idx) );
            if strcmp(DynamicsModel,'THode')
                Crtl_Aug = A_inv*[zeros(3,1); Crtl(jj).u(:,i)];
                Dim_Crtl(jj).u(:,i) = Crtl_Aug(4:6);
            elseif strcmp(DynamicsModel,'CWode')
                Dim_Crtl(jj).u(:,i) = Tc*Crtl(jj).u(:,i);
            elseif strcmp(DynamicsModel,'NLode')
                Dim_Crtl(jj).u(:,i) = Rrho*Crtl(jj).u(:,i);
            end
        end
    end
    
    if exist('u11','file') == 3
        delete u11.c u11.mexmaci64
    end
    if exist('u12','file') == 3
        delete u12.c u12.mexmaci64
    end
    if exist('u13','file') == 3
        delete u13.c u13.mexmaci64
    end
    if exist('u21','file') == 3
        delete u21.c u21.mexmaci64
    end
    if exist('u22','file') == 3
        delete u22.c u22.mexmaci64
    end
    if exist('u23','file') == 3
        delete u23.c u23.mexmaci64
    end
    
    u11 = casadi.interpolant('U','bspline',{xnew}, Crtl(1).u(1,:));
    u12 = casadi.interpolant('U','bspline',{xnew}, Crtl(1).u(2,:));
    u13 = casadi.interpolant('U','bspline',{xnew}, Crtl(1).u(3,:));
    u21 = casadi.interpolant('U','bspline',{xnew}, Crtl(2).u(1,:));
    u22 = casadi.interpolant('U','bspline',{xnew}, Crtl(2).u(2,:));
    u23 = casadi.interpolant('U','bspline',{xnew}, Crtl(2).u(3,:));
    
    u11.generate('u11.c',CasADiopts);
    u12.generate('u12.c',CasADiopts);
    u13.generate('u13.c',CasADiopts);
    u21.generate('u21.c',CasADiopts);
    u22.generate('u22.c',CasADiopts);
    u23.generate('u23.c',CasADiopts);
    mex u11.c -largeArrayDims
    mex u12.c -largeArrayDims
    mex u13.c -largeArrayDims
    mex u21.c -largeArrayDims
    mex u22.c -largeArrayDims
    mex u23.c -largeArrayDims
    clear u11 u12 u13 u21 u22 u23
end
% -------------------------------------------------------------------------

%% ------------------ Integrate the system's trajectory -------------------
for ll = 1
    disp('Integrating the resulting trajectory...');
    if strcmp(DynamicsModel,'THode')
        [~,X_ode45] = ode113(@(t,x)TH_ode(t,x),xnew,X0,options);
    elseif strcmp(DynamicsModel,'CWode')
        [~,X_ode45] = ode113(@(t,x)CW_ode(t,x),xnew,X0,options);
    elseif strcmp(DynamicsModel,'NLode')
        [~,X_ode45_1] = ode113(@(t,x)NL_ode(t,x,1),xnew,X0(1:6),options);
        [~,X_ode45_2] = ode113(@(t,x)NL_ode(t,x,2),xnew,X0(7:12),options);
        X_ode45 = [X_ode45_1 X_ode45_2];
    end
    
    if strcmp(DynamicsModel,'THode')
        for i = 1:length(xnew)
            f      = xnew(i);
            h      = sqrt(mu_Earth*a_chief*(1-e_chief^2));
            K_f    = 1+e_chief*cos(f);
            conts1 = h^(3/2)/(mu_Earth*K_f);
            conts2 = mu_Earth*e_chief*sin(f)/h^(3/2);
            conts3 = 1/conts1;
            A_inv  = [conts1*eye(3)       zeros(3);
                conts2*eye(3)  conts3*eye(3)];
            X_ode45_LVLH(i,:) = (blkdiag(A_inv,A_inv)*X_ode45(i,:).').';
            Pint(:,i)         = blkdiag(A_inv,A_inv)*Pint(:,i);
        end
    elseif strcmp(DynamicsModel,'CWode')
        X_ode45_LVLH = Tc*X_ode45;
        Pint         = Tc*Pint;
    elseif strcmp(DynamicsModel,'NLode')
        X_ode45_LVLH = [(blkdiag(Rrho,Rrho)*X_ode45_1.').' (blkdiag(Rrho,Rrho)*X_ode45_2.').'];
        Pint         = [blkdiag(Rrho,Rrho)*Pint(1:6,:); blkdiag(Rrho,Rrho)*Pint(7:12,:)];
    end
    
    disp('Done!!!!!');
    
    % Solution.ConvergedTrajectory = X_ode45_LVLH;
    % Solution.Agent1RequiredControl = Crtl(1).u;
    % Solution.Agent2RequiredControl = Crtl(2).u;
    % Solution.CostFunction = cost;
    % save('MinimumFuel_Geometric_Planner.mat', 'Solution')
end

%% ------------------------------------------------------------------------
% ######################## Plotting Results ###############################
% -------------------------------------------------------------------------
%% ------------------- Plotting the extracted control ---------------------
for ll = 1
    fh2 = figure;
    Label  = {'Agent 1', 'Agent 2', 'Agent 3', 'Agent 4'};
    for i = 1:Num_Agent
        
        if strcmp(DynamicsModel,'THode')
            Time = (Tc/(2*pi))*xnew/3600;
        else
            Time = Tc*xnew/3600;
        end
        set(gca,'XTick',0:pi/2:2*pi)
        set(gca,'XTickLabel',{'0','$\pi$/2','$\pi$','3$\pi$/2','2$\pi$'})
        subplot(2,2,i)
        plot(Time,Dim_Crtl(i).u(1,:),'r','LineWidth',2.5);
        hold on
        plot(Time,Dim_Crtl(i).u(2,:),'b','LineWidth',2.5);
        plot(Time,Dim_Crtl(i).u(3,:),'g','LineWidth',2.5);
        ylabel({'Required Control', '(LVLH)'})
        xlabel('time (hr)')
        xlim([Time(1) Time(end)])
        grid on;
        title(Label(i))
        
        
        subplot(2,2,i+2)
        plt1 = plot(xnew,Crtl(i).u(1,:),'r','LineWidth',2.5);
        hold on
        plt2 = plot(xnew,Crtl(i).u(2,:),'b','LineWidth',2.5);
        plt3 = plot(xnew,Crtl(i).u(3,:),'g','LineWidth',2.5);
        ylabel({'Nondimensionalized', 'control'})
        if strcmp(DynamicsModel,'THode')
            xlabel('f [rad]')
        else
            xlabel('$\tau [-]$')
        end
        xlim([xnew(1) xnew(end)])
        grid on;
        title(Label(i))
    end
    hL = legend([plt1, plt2, plt3],...
        {'u1','u2','u3'},'AutoUpdate','off');
    hL.FontSize = 20;
    newPosition = [0.22 0.14 0.1 0.1];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1  );
    
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do not show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 800 700];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 5580 4320]/res;
    % save figure
    % print(fh2,'./Figures/RequiredControl','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% --------- Plotting the action integral across each iteratons ----------
for ll = 1
    fh2 = figure;
    loglog(t,cost,'Color',c5,'LineWidth',4);
    title('Action integral across iteration')
    xlabel('Homotopy iteration')
    ylabel('Length of the actuated curve')
    grid on;
    set(gca,'FontSize',30)
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 700 600];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 5580 4320]/res;
    % save figure
    % print(fh2,'./Figures/ActionIntegral','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% ---- 3D plot of the converged solution (Nondimensionalized units) ------
for ll = 1
    % cMat1 = [c5;c6];
    % cMat2 = [c3;c1];
    % cMat0 = [c7;c4];
    %
    % fh2 = figure;
    % hold on
    % plt0  = zeros(Num_Agent);
    % plt1   = zeros(Num_Agent);
    % plt2  = zeros(Num_Agent);
    % for jj = 1:Num_Agent
    %     idx = 1+mm*(jj-1):mm*jj;
    %
    %     plt0(:,jj) = plot3(X_ode45(:,idx(1)),X_ode45(:,idx(2)),X_ode45(:,idx(3)),...
    %         '--','Color',cMat0(jj,:),'LineWidth',5);
    %
    %     plt1(:,jj) = plot3(X_rel_0(1,idx(1)),X_rel_0(1,idx(2)),X_rel_0(1,idx(3)),'*',...
    %         'LineWidth',3,...
    %         'MarkerEdgeColor',cMat1(jj,:),...
    %         'MarkerFaceColor',cMat1(jj,:)',...
    %         'MarkerSize',8);
    %
    %
    %     plt2(:,jj) = plot3(X_rel_1(1,idx(1)),X_rel_1(1,idx(2)),X_rel_1(1,idx(3)),'o',...
    %         'LineWidth',3,...
    %         'MarkerEdgeColor',cMat2(jj,:),...
    %         'MarkerFaceColor',cMat2(jj,:),...
    %         'MarkerSize',8);
    %
    %     h5 = plot3(0,0,0,'bo',...
    %         'LineWidth',2,...
    %         'MarkerEdgeColor','k',...
    %         'MarkerFaceColor','k',...
    %         'MarkerSize',15);
    %     plot3(X_rel_0(:,idx(1)),X_rel_0(:,idx(2)),X_rel_0(:,idx(3)),'Color',cMat1(jj,:),'LineWidth',5);
    %     plot3(X_rel_1(:,idx(1)),X_rel_1(:,idx(2)),X_rel_1(:,idx(3)),'Color',cMat2(jj,:),'LineWidth',5);
    %
    %     grid on
    %     xlabel('$\bar{x}$')
    %     ylabel('$\bar{y}$')
    %     zlabel('$\bar{z}$')
    %     title('3D Relative Trajectory (Nondimensionalized units)')
    % end
    % % r = 1;
    % % [x_sphere, y_sphere, z_sphere] = sphere(100);
    % % h = surf(r*x_sphere, r*y_sphere, r*z_sphere);
    % % set(h,'FaceAlpha', 0.09, 'EdgeColor', 'b', 'EdgeAlpha',.05, 'FaceColor', 'b');
    % % axis equal
    % set(gca,'FontSize',30)
    % view(68,10)
    % hL = legend([h5, plt1(end,1), plt1(end,2), plt2(end,1), plt2(end,2), plt0(end,1), plt0(end,2)],...
    %     {'Chief Location','Agent1 Init Condition','Agent2 Init Condition',...
    %     'Agent1 End Condition','Agent2 End Condition','Agent1 Transfer Traj',...
    %     'Agent2 Transfer Traj'},'AutoUpdate','off');
    % % Programatically move the Legend
    % newPosition = [0.21 0.72 0.1 0.1];
    % newUnits = 'normalized';
    % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    %
    % % set all units inside figure to normalized so that everything is scaling accordingly
    % set(findall(fh2,'Units','pixels'),'Units','normalized');
    % % do show figure on screen
    % set(fh2, 'visible', 'on')
    % % set figure units to pixels & adjust figure size
    % fh2.Units = 'pixels';
    % fh2.OuterPosition = [10 10 900 800];
    % % define resolution figure to be saved in dpi
    % res = 500;
    % % recalculate figure size to be saved
    % set(fh2,'PaperPositionMode','manual')
    % fh2.PaperUnits = 'inches';
    % fh2.PaperPosition = [0 0 5580 4320]/res;
    % % save figure
    % % print(fh2,'./Figures/NondimensionalizedConvergeSolution','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% --- 3D plot of the converged solution (in the configuration space) -----
for ll = 1
    cMat1 = [c5;c6];
    cMat2 = [c3;c1];
    cMat0 = [c7;c4];
    fh2 = figure;
    hold on
    plt0  = zeros(Num_Agent);
    plt1   = zeros(Num_Agent);
    plt2  = zeros(Num_Agent);
    
    for jj = 1:Num_Agent
        idx = 1+mm*(jj-1):mm*jj;
        plt0(:,jj) = plot3(X_ode45_LVLH(:,idx(1)),X_ode45_LVLH(:,idx(2)),X_ode45_LVLH(:,idx(3)),...
            '-.','Color',cMat0(jj,:),'LineWidth',4);
        
        plt1(:,jj) = plot3(X_THpos01(1,idx(1)),X_THpos01(1,idx(2)),X_THpos01(1,idx(3)),'*',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat1(jj,:),...
            'MarkerFaceColor',cMat1(jj,:)',...
            'MarkerSize',8);
        
        plt2(:,jj) = plot3(X_THpos11(1,idx(1)),X_THpos11(1,idx(2)),X_THpos11(1,idx(3)),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat2(jj,:),...
            'MarkerFaceColor',cMat2(jj,:),...
            'MarkerSize',8);
        
        h5 = plot3(0,0,0,'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','k',...
            'MarkerSize',15);
        
        plot3(X_THpos11(:,idx(1)),X_THpos11(:,idx(2)),X_THpos11(:,idx(3)),'Color',cMat2(jj,:),'LineWidth',4);
        plot3(X_THpos01(:,idx(1)),X_THpos01(:,idx(2)),X_THpos01(:,idx(3)),'Color',cMat1(jj,:),'LineWidth',4);
        
        grid on
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('3D Relative Trajectory')
    end
    
    set(gca,'FontSize',30)
    view(68,10)
    % hL = legend([h5, plt1(end,1), plt1(end,2), plt2(end,1), plt2(end,2), plt0(end,1), plt0(end,2)],...
    %     {'Chief Location','Agent1 Init Condition','Agent2 Init Condition',...
    %     'Agent1 End Condition','Agent2 End Condition','Agent1 Transfer Traj',...
    %     'Agent2 Transfer Traj'},'AutoUpdate','off');
    % % Programatically move the Legend
    % % saveLegendToImage(fh2, hL, './Figures/3DLegend')
    % newPosition = [0.22 0.72 0.1 0.1]; %[0.21 0.72 0.1 0.1];
    % newUnits = 'normalized';
    % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 900 800];
    % define resolution figure to be saved in dpi
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 5580 4320]/res;
    % save figure
    % print(fh2,'./Figures/ConvergeSolution','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% ------ Animation of the homotopy trasformation during convergence ------
for k = 1
    % cMat1  = [c5;c6];
    % cMat2 = [c3;c1];
    % cMat0 = [c7;c4];
    %
    % fh2 = figure;
    % view(68,10)
    % hold on
    % plt0  = zeros(Num_Agent);
    % plt1   = zeros(Num_Agent);
    % plt2  = zeros(Num_Agent);
    % for jj = 1%:Num_Agent
    %     idx = 1+mm*(jj-1):mm*jj;
    %
    %     plt0(:,jj) = plot3(X_ode45(:,idx(1)),X_ode45(:,idx(2)),X_ode45(:,idx(3)),...
    %        'Color',cMat0(jj,:),'LineWidth',4);
    %
    %     plt1(:,jj) = plot3(X_rel_0(1,idx(1)),X_rel_0(1,idx(2)),X_rel_0(1,idx(3)),'*',...
    %         'LineWidth',3,...
    %         'MarkerEdgeColor',cMat1(jj,:),...
    %         'MarkerFaceColor',cMat1(jj,:)',...
    %         'MarkerSize',8);
    %
    %
    %     plt2(:,jj) = plot3(X_rel_1(1,idx(1)),X_rel_1(1,idx(2)),X_rel_1(1,idx(3)),'o',...
    %         'LineWidth',3,...
    %         'MarkerEdgeColor',cMat2(jj,:),...
    %         'MarkerFaceColor',cMat2(jj,:),...
    %         'MarkerSize',8);
    %
    %     h5 = plot3(0,0,0,'bo',...
    %         'LineWidth',2,...
    %         'MarkerEdgeColor','k',...
    %         'MarkerFaceColor','k',...
    %         'MarkerSize',15);
    %     h3 = plot3(X_rel_0(:,idx(1)),X_rel_0(:,idx(2)),X_rel_0(:,idx(3)),'-.','Color',cMat1(jj,:),'LineWidth',4);
    %     h4 = plot3(X_rel_1(:,idx(1)),X_rel_1(:,idx(2)),X_rel_1(:,idx(3)),'-.','Color',cMat2(jj,:),'LineWidth',4);
    % end
    %
    % % set all units inside figure to normalized so that everything is scaling accordingly
    % set(findall(fh2,'Units','pixels'),'Units','normalized');
    % % do show figure on screen
    % set(fh2, 'visible', 'on')
    % % set figure units to pixels & adjust figure size
    % fh2.Units = 'pixels';
    % fh2.OuterPosition = [10 10 900 800];
    % % define resolution figure to be saved in dpi
    % res = 500;
    % % recalculate figure size to be saved
    % set(fh2,'PaperPositionMode','manual')
    % fh2.PaperUnits = 'inches';
    % fh2.PaperPosition = [0 0 5580 4320]/res;
    %
    %
    % tpointsCorase = 250;
    % xpointsCorase = 250;
    % Tmesh = [0, logspace(-4,log10(tmax),tpointsCorase-1)]; % discretization of time
    % Xmesh = linspace(0,T,xpointsCorase); % discretization of the curve
    % [X,Y] = meshgrid(Tmesh,Xmesh);
    %
    % Agent1_u1 = interp2(t,x,sol(:,:,1),X,Y);
    % Agent1_u2 = interp2(t,x,sol(:,:,2),X,Y);
    % Agent1_u3 = interp2(t,x,sol(:,:,3),X,Y);
    %
    % Xrel = Rrho*[Agent1_u1(1,:);Agent1_u2(1,:);Agent1_u3(1,:)];
    % X = Xrel(1,:);
    % Y = Xrel(2,:);
    % Z = Xrel(3,:);
    % h2 = plot3(X,Y,Z,'Color',c9,'LineWidth',4);
    % hold on
    % h1 = plot3(X,Y,Z,'r','LineWidth',4);
    % grid on;
    % title('Nonlinear Geomeric Planning - 3D configuration');
    % xlabel('$x$');
    % ylabel('$y$');
    % zlabel('$z$');
    % hL = legend([plt0(end,1), h2, h1,h3, h4, h5],...
    %     {'Converged Solution','Initial Guess','Homotopy Iteration','Departure Orbit','Arrival Orbit','Chief'},'AutoUpdate','off');
    % hL.FontSize = 27;
    % newPosition = [.3 .775 0.01 0.01];
    % newUnits = 'normalized';
    % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    % set(gca,'nextplot','replacechildren','visible','off')
    % f  = getframe;
    % [im,map] = rgb2ind(f.cdata,1000,'nodither'); % 65536
    % im(1,1,1,tpointsCorase) = 0;
    % for i = 1:tpointsCorase
    %     Xrel = Rrho*[Agent1_u1(i,:);Agent1_u2(i,:);Agent1_u3(i,:)];
    %     X = Xrel(1,:); h1.XDataSource='X';
    %     Y = Xrel(2,:); h1.YDataSource='Y';
    %     Z = Xrel(3,:); h1.ZDataSource='Z';
    %     refreshdata(h1,'caller');
    %     drawnow;
    %     f = getframe;
    %     im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
    % end
    % imwrite(im,map,'./Figures/AGHF_THode.gif','DelayTime',0,'LoopCount',inf) %g443800
end
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% ######################## Required Functions #############################
% -------------------------------------------------------------------------
%% --------------------------- PDE for pdepe ------------------------------
function [c,f,s] = mypdexpde(x,t,u,DuDx,k,N) % Define PDE; right-hand-side of AGHF

XChief = full([Chief_x(x);  Chief_y(x);  Chief_z(x);...
               Chief_dx(x); Chief_dy(x); Chief_dz(x)]);
LL = full(L(u,DuDx,k,x));
CasADiresult = full(dLdx(u,DuDx,k,x,XChief,LL))';
CasADir = [CasADiresult(1:12), CasADiresult(13:24)];

% dL/dx
pLx = CasADir(:,1); % EL(:,1); %
% dL/d(dot_x)
pLxd = CasADir(:,2); % EL(:,2); % 

f = pLxd;
s = -pLx;
c = ones(N,1);

end
% -------------------------------------------------------------------------

%% --------------------- Initial condition for pdepe ----------------------
function u0 = mypdexic(x, X0, Xf, T)    % Initial condition for AGHF
%%  straight line connecting X0 and Xf
% u0=X0+(Xf-X0)*(x/T);
% % add some sinusoidal deviation to x1 initial guess (can be discarded)
% u0(1)=X0(1)+(Xf(1)-X0(1))*(x/T) + 0.01*sin(2*pi*x/T);

%% Sinusoidal initial condition
freq1 = pi/2 + 2*pi;
freq2 = pi + 2*pi;

u0 = X0*cos(freq1*x/T) ...
    + Xf*sin(freq1*x/T) ...
    - (X0+Xf)*sin(freq2*x/T);

end
% -------------------------------------------------------------------------

%% --------------------- Boundary condition for pdepe ---------------------
function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t, X0, Xf, N) % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end
% -------------------------------------------------------------------------

function [Xdot] = NonDimentionalized_NLode(t,X)
global mu_Earth Rrho Tc 
format long
X_chief = full([Chief_x(t);  Chief_y(t);  Chief_z(t);...
                Chief_dx(t); Chief_dy(t); Chief_dz(t)]);
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
rt = X_chief(1:3);  vt = X_chief(4:6); rt_norm = (rt.'*rt).^(1/2);
rho_bar = X(1:3,:); drho_bar = X(4:6,:);

% Redimensionalizing the variables in the ODE
rho = Rrho*rho_bar; rho_prime = Rrho*drho_bar;
rc  = [rt_norm; 0; 0] + rho; rc_norm =  (rc.'*rc).^(1/2);

% Calculating the respective accelerations
TN = DCM(rt,vt); % DCM's between target's RIC frame and ECI frame

% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(rt)*vt;
Omega = TN*( h_vec/rt_norm^2);
Omegadot = TN*( -2*(rt.'*vt)*h_vec/rt_norm^4 );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
    - tilde(Omega)*(tilde(Omega)*rho) ...
    - tilde(Omegadot)*rho;

% Nondimetional relative ODE of the deputy
Xdot = [drho_bar; (del_ag + Coriolis_Effect)];
Xdot = (blkdiag(Rrho,Rrho)\Xdot)*Tc;

end
function [Xdot] = NonDimentionalized_NLode0(t,X)
global mu_Earth Rrho Tc 
format long
X_chief = full([Chief_x0(t);  Chief_y0(t);  Chief_z0(t);...
                Chief_dx0(t); Chief_dy0(t); Chief_dz0(t)]);
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
rt = X_chief(1:3);  vt = X_chief(4:6); rt_norm = (rt.'*rt).^(1/2);
rho_bar = X(1:3,:); drho_bar = X(4:6,:);

% Redimensionalizing the variables in the ODE
rho = Rrho*rho_bar; rho_prime = Rrho*drho_bar;
rc  = [rt_norm; 0; 0] + rho; rc_norm =  (rc.'*rc).^(1/2);

% Calculating the respective accelerations
TN = DCM(rt,vt); % DCM's between target's RIC frame and ECI frame

% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(rt)*vt;
Omega = TN*( h_vec/rt_norm^2);
Omegadot = TN*( -2*(rt.'*vt)*h_vec/rt_norm^4 );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
    - tilde(Omega)*(tilde(Omega)*rho) ...
    - tilde(Omegadot)*rho;

% Nondimetional relative ODE of the deputy
Xdot = [drho_bar; (del_ag + Coriolis_Effect)];
Xdot = (blkdiag(Rrho,Rrho)\Xdot)*Tc;

end
function [Xdot] = NonDimentionalized_NLodef(t,X)
global mu_Earth Rrho Tc 
format long
X_chief = full([Chief_xf(t);  Chief_yf(t);  Chief_zf(t);...
                Chief_dxf(t); Chief_dyf(t); Chief_dzf(t)]);
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
rt = X_chief(1:3);  vt = X_chief(4:6); rt_norm = (rt.'*rt).^(1/2);
rho_bar = X(1:3,:); drho_bar = X(4:6,:);

% Redimensionalizing the variables in the ODE
rho = Rrho*rho_bar; rho_prime = Rrho*drho_bar;
rc  = [rt_norm; 0; 0] + rho; rc_norm =  (rc.'*rc).^(1/2);

% Calculating the respective accelerations
TN = DCM(rt,vt); % DCM's between target's RIC frame and ECI frame

% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(rt)*vt;
Omega = TN*( h_vec/rt_norm^2);
Omegadot = TN*( -2*(rt.'*vt)*h_vec/rt_norm^4 );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
    - tilde(Omega)*(tilde(Omega)*rho) ...
    - tilde(Omegadot)*rho;

% Nondimetional relative ODE of the deputy
Xdot = [drho_bar; (del_ag + Coriolis_Effect)];
Xdot = (blkdiag(Rrho,Rrho)\Xdot)*Tc;

end
% -------------------------------------------------------------------------

function [Xdot] = NonDimentionalized_CWode(t,X)
global mu_dot Tc

% Computing the CW control and dynamic matrices
A  = [zeros(3) eye(3);
    3*mu_dot^2 0 0 0 2*mu_dot 0;
    0 0 0 -2*mu_dot 0 0;
    0 0 -mu_dot^2 0 0 0]*Tc;

% Reshaping the input into an 6 by n matrix
[l,~] = size(X);
X = reshape(X,6,l/6);
Xdot = A*X;
Xdot = Xdot(:);
end
% -------------------------------------------------------------------------

function [Xdot] = NonDimentionalized_THode(t,X)
global e_chief
k = 1+e_chief*cos(t);
A = [zeros(3) eye(3);
    3/k  0  0   0   2  0;
    0    0  0   -2  0  0;
    0    0  -1  0   0  0];
[l,~] = size(X);
X = reshape(X,6,l/6);
Xdot = A*X;
Xdot = Xdot(:);
end
% -------------------------------------------------------------------------

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!   CasADi Functions   !!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function [Xdot] = NonDimentionalized_NLodeCasADi(t,X,XChief)
global mu_Earth Rrho Tc
format long 
tilde   = @(v) [0    -v(3)  v(2);
                v(3)   0   -v(1);
               -v(2)  v(1)    0];
% Collecting relevant quantities
rt = XChief(1:3); vt = XChief(4:6); rt_norm = (rt.'*rt).^(1/2);
rho_bar = X(1:3,:); drho_bar = X(4:6,:);

% Redimensionalizing the variables in the ODE
rho = Rrho*rho_bar; rho_prime = Rrho*drho_bar;
rc = [rt_norm; 0; 0] + rho; rc_norm =  (rc.'*rc).^(1/2);

% Calculating the respective accelerations
TN = DCM(rt,vt); % DCM's between target's RIC frame and ECI frame

% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(rt)*vt;
Omega = TN*( h_vec/rt_norm^2);
Omegadot = TN*( -2*(rt.'*vt)*h_vec/rt_norm^4 );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
    - tilde(Omega)*(tilde(Omega)*rho) ...
    - tilde(Omegadot)*rho;

% Nondimetional relative ODE of the deputy
Xdot = [drho_bar; (del_ag + Coriolis_Effect)];
Xdot = (blkdiag(Rrho,Rrho)\Xdot)*Tc;

end
% -------------------------------------------------------------------------

function [Xdot] = NonDimentionalized_CWodeCasADi(t,X,XChief)
global mu_dot Tc

% Computing the CW control and dynamic matrices
A  = [zeros(3) eye(3);
    3*mu_dot^2 0 0 0 2*mu_dot 0;
    0 0 0 -2*mu_dot 0 0;
    0 0 -mu_dot^2 0 0 0]*Tc;

% Reshaping the input into an 6 by n matrix
[l,~] = size(X);
X = reshape(X,6,l/6);
Xdot = A*X;
Xdot = Xdot(:);
end
% -------------------------------------------------------------------------

function [Xdot] = NonDimentionalized_THodeCasADi(t,X,XChief)
global e_chief
k = 1+e_chief*cos(t);
A = [zeros(3) eye(3);
    3/k  0  0   0   2  0;
    0    0  0   -2  0  0;
    0    0  -1  0   0  0];
[l,~] = size(X);
X = reshape(X,6,l/6);
Xdot = A*X;
Xdot = Xdot(:);
end
% -------------------------------------------------------------------------



