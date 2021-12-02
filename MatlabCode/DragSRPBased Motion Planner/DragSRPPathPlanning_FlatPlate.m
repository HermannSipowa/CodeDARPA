%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hermann Kaptui Sipowa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
clear
close all
clc
start_up
format long e
global mu_Earth Target Chasser ...
       Tc Num_Agent ae JD AU % iteration InitalCondition
ae        = 6378.136;
mu_Earth  = 3.986004415E5;
JD        = 2456296.25;
AU        = 149597870.7;
Chasser   = Spacecraft3D([20, 20, 1, .25, 2, 2, 0.2, 0.5, 3, 0, 0]); % Deputy characteristic
Target    = Spacecraft3D([0, 0, 2, .25, 1.5, .5, 0.2, 0.5, 1.5, 0, 0]); % Target characteristic
Num_Agent = 1; % Number of agents in the system
SimTime   = 0.25;
N   = 13;
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

%% Calculating the boundaries conditions
%=======================================%
% Specifying the chief's orbit
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
hp           = 1e4; % Peri-apsis altitude
e_chief      = 0.1; % Eccentricity
a_chief      = (ae+hp)/(1-e_chief); % Semi-major axis in Km 
inc_chief    = 50;  % Inclination in deg0
BigOmg_chief = 10;  % RAAN in deg
LitOmg_chief = 10;  % AOP in deg
M_chief1     = 0;   % Initial mean anomaly

% a_chief      = 1e4; % Semi-major axis in Km
% e_chief      = 0.5; % Eccentricity
% inc_chief    = 50;  % Inclination in deg0
% BigOmg_chief = 10;  % RAAN in deg
% LitOmg_chief = 10;  % AOP in deg
% M_chief1     = 0;   % Initial mean anomaly

%% Integrating the chief's trajectory
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
for i = 1
    
    Period  = (2*pi)*sqrt(a_chief^3/mu_Earth);
    Tc      = 4*Period;
    IntTime = SimTime*Tc;
    Tsize   = 1e4;
    tspan   = linspace(0,IntTime,Tsize);
    options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-30);
    mu_dot  = sqrt(mu_Earth/a_chief^3); % Chief's mean motion
    Tnorm   = tspan/Tc;
    Ttraj1  = linspace(0,Period,Tsize); TT1 = Ttraj1/Tc;
    Ttraj2  = IntTime+linspace(0,Period,Tsize); TT2 = Ttraj2/Tc;
    
    inc_chief = deg2rad(inc_chief); BigOmg_chief = deg2rad(BigOmg_chief);
    LitOmg_chief = deg2rad(LitOmg_chief); M_chief1 = deg2rad(M_chief1);
    COE1 = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief1];
    [Position_target1,Velocity_target1] = COEstoRV(COE1,mu_Earth);
    X0_Chief1  = [Position_target1; Velocity_target1];
    
    [~, Xnom]  = ode113(@(t,X)DragSRPChiefMotionODE(t,X,Target),tspan,X0_Chief1,options);
    [~, Xnom0] = ode113(@(t,X)DragSRPChiefMotionODE(t,X,Target),Ttraj1,Xnom(1,:)',options);
    [~, Xnomf] = ode113(@(t,X)DragSRPChiefMotionODE(t,X,Target),Ttraj2,Xnom(end,:)',options);
end

%% ------- Create MX files to interpolate the chief's states (Xnom) -------
for i = 1
    addpath('../casadi-osx-matlabR2015a-v3.5.5')
    import casadi.*
    CasADiopts = struct('main', true, 'mex', true);
    tt = MX.sym('t');
    %=====================================================================%
    if exist('ChiefState','file') == 3
        delete ChiefState.c ChiefState.mexmaci64
    end
    Chief_x  = casadi.interpolant('Chief_x','bspline',{Tnorm},Xnom(:,1)');
    Chief_y  = casadi.interpolant('Chief_y','bspline',{Tnorm},Xnom(:,2)');
    Chief_z  = casadi.interpolant('Chief_z','bspline',{Tnorm},Xnom(:,3)');
    Chief_dx = casadi.interpolant('Chief_dx','bspline',{Tnorm},Xnom(:,4)');
    Chief_dy = casadi.interpolant('Chief_dy','bspline',{Tnorm},Xnom(:,5)');
    Chief_dz = casadi.interpolant('Chief_dz','bspline',{Tnorm},Xnom(:,6)');
    Output   = [Chief_x(tt);Chief_y(tt);Chief_z(tt);...
                Chief_dx(tt);Chief_dy(tt);Chief_dz(tt)];
    ChiefState = Function('ChiefState',{tt},...
                 {Output},{'time'},{'Output'});
    ChiefState.generate('ChiefState.c',CasADiopts);
    mex ChiefState.c -largeArrayDims
    
    
    %=====================================================================%
    if exist('ChiefState0','file') == 3
        delete ChiefState0.c ChiefState0.mexmaci64
    end
    Chief_x0  = casadi.interpolant('Chief_x0','bspline',{TT1},Xnom0(:,1)');
    Chief_y0  = casadi.interpolant('Chief_y0','bspline',{TT1},Xnom0(:,2)');
    Chief_z0  = casadi.interpolant('Chief_z0','bspline',{TT1},Xnom0(:,3)');
    Chief_dx0 = casadi.interpolant('Chief_dx0','bspline',{TT1},Xnom0(:,4)');
    Chief_dy0 = casadi.interpolant('Chief_dy0','bspline',{TT1},Xnom0(:,5)');
    Chief_dz0 = casadi.interpolant('Chief_dz0','bspline',{TT1},Xnom0(:,6)');
    Output0   = [Chief_x0(tt);Chief_y0(tt);Chief_z0(tt);...
                 Chief_dx0(tt);Chief_dy0(tt);Chief_dz0(tt)];
    ChiefState0 = Function('ChiefState0',{tt},...
                  {Output0},{'time'},{'Output0'});
    ChiefState0.generate('ChiefState0.c',CasADiopts);
    mex ChiefState0.c -largeArrayDims
    
    
    %=====================================================================%
    if exist('ChiefStatef','file') == 3
        delete ChiefStatef.c ChiefStatef.mexmaci64
    end
    Chief_xf  = casadi.interpolant('Chief_xf','bspline',{TT2},Xnomf(:,1)');
    Chief_yf  = casadi.interpolant('Chief_yf','bspline',{TT2},Xnomf(:,2)');
    Chief_zf  = casadi.interpolant('Chief_zf','bspline',{TT2},Xnomf(:,3)');
    Chief_dxf = casadi.interpolant('Chief_dxf','bspline',{TT2},Xnomf(:,4)');
    Chief_dyf = casadi.interpolant('Chief_dyf','bspline',{TT2},Xnomf(:,5)');
    Chief_dzf = casadi.interpolant('Chief_dzf','bspline',{TT2},Xnomf(:,6)');
    Outputf   = [Chief_xf(tt);Chief_yf(tt);Chief_zf(tt);...
                 Chief_dxf(tt);Chief_dyf(tt);Chief_dzf(tt)];
    ChiefStatef = Function('ChiefStatef',{tt},...
                  {Outputf},{'time'},{'Outputf'});
    ChiefStatef.generate('ChiefStatef.c',CasADiopts);
    mex ChiefStatef.c -largeArrayDims
    
    %=====================================================================%
    clear Chief_x Chief_y Chief_z Chief_dx Chief_dy Chief_dz Output ChiefState ...
          Chief_x0 Chief_y0 Chief_z0 Chief_dx0 Chief_dy0 Chief_dz0 Output0 ChiefState0 ...
          Chief_xf Chief_yf Chief_zf Chief_dxf Chief_dyf Chief_dzf Outputf ChiefStatef
end
% -------------------------------------------------------------------------

%% --------------- Setting the deputy initial conditions ------------------
for ll = 1
    % % Deputy initial conditions
    % %~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % r0 = Xnom0(1,1:3)'; v0 = Xnom0(1,4:6)';
    % COE0 = RVtoCOEs(r0,v0,mu_Earth); f1 = COE0(6);
    % q1_chief = e_chief*cos(LitOmg_chief); q2_chief = e_chief*sin(LitOmg_chief);
    % theta_chief1 = f1+LitOmg_chief;
    % OE_chief1 = [a_chief, theta_chief1, inc_chief, q1_chief, q2_chief, BigOmg_chief];
    % AMap1 = ForwardMapping(OE_chief1, mu_Earth); % Linear mapping matrix
    %
    % dela = 0;
    % dele = 25/(5*a_chief);
    % deli = -50/(3*a_chief);
    % delLitOmg = -2*pi*1E-6;
    % delBigOmg = 0;
    % delM = pi*1E-4;
    % delq1_1 = dele*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    % delq2_1 = dele*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    % deltheta = delLitOmg + delM; % Relative true latitude in rad
    % delCOE_1 = [dela, deltheta, deli, delq1_1, delq2_1, delBigOmg]';
    % X01 = AMap1*delCOE_1;
    %
    % % Deputy final conditions
    % %~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % af_Deyty      = a_chief+25; % Semi-major axis in Km (600 Km altitude)
    % ef_Deyty      = 0.15;       % Eccentricity
    % incf_Deyty    = 50.1;       % Inclination in deg0
    % BigOmgf_Deyty = 10.0;       % RAAN in deg
    % LitOmgf_Deyty = 10.0;          % AOP in deg
    % Mf_Deyty      = 0;          % Initial mean anomaly
    % incf_Deyty = deg2rad(incf_Deyty); BigOmgf_Deyty = deg2rad(BigOmgf_Deyty);
    % LitOmgf_Deyty = deg2rad(LitOmgf_Deyty); Mf_Deyty = deg2rad(Mf_Deyty);
    % COEDeputyfinal = [af_Deyty,ef_Deyty,incf_Deyty,BigOmgf_Deyty,LitOmgf_Deyty,Mf_Deyty];
    % [FinalDeputyPosition,DeputyFinalVelocity] = COEstoRV(COEDeputyfinal,mu_Earth);
    %
    % r_target = Xnom(end,1:3)'; v_target = Xnom(end,4:6)';
    % TN = DCM(r_target,v_target); NOmega_target = cross(r_target,v_target)/norm(r_target)^2;
    % NR_rel = FinalDeputyPosition - r_target; NV_rel = DeputyFinalVelocity - v_target;
    % rho    = TN*NR_rel; rhodot = TN*(NV_rel - cross(NOmega_target,NR_rel));
    % Xf1    = [rho; rhodot];
    %
    %
    % % q_t0 = load('q_t0.mat');  qr_chasser00 = q_t0.qr_chasser00;
    % % q_tf = load('q_tf.mat');  qr_chasserf0 = q_tf.qr_chasserf0;
    % qr_chasser00 = [6; 4; 9; 1]; % randi([1 9],[4 1]); %
    % qr_chasserf0 = [7; 5; 1; 8]; % randi([1 9],[4 1]); %
    %
    % qr_chasser0 = qr_chasser00/norm(qr_chasser00);
    % Omega_chaser0 = 5e-2*[deg2rad(-.9) deg2rad(.8) -deg2rad(.7)]'; % Angular velocity in inertial frame
    % CN = Quaternion_to_DCM(qr_chasser0); Omega_chaser_B0 = CN*Omega_chaser0;
    % X0 = [qr_chasser0; Omega_chaser0; X01];
    %
    % qr_chasserf = qr_chasserf0/norm(qr_chasserf0);
    % Omega_chaserf = 3e-3*[deg2rad(-.3) deg2rad(.3) -deg2rad(.1)]';
    % CN = Quaternion_to_DCM(qr_chasserf); Omega_chaser_Bf = CN*Omega_chaserf;
    % Xf = [qr_chasserf; Omega_chaserf; Xf1];
end
for ll = 1
    r0 = Xnom0(1,1:3)'; v0 = Xnom0(1,4:6)';
    COE0 = RVtoCOEs(r0,v0,mu_Earth); f1 = COE0(6);
    q1_chief = e_chief*cos(LitOmg_chief); q2_chief = e_chief*sin(LitOmg_chief);
    theta_chief1 = f1+LitOmg_chief;
    OE_chief1 = [a_chief, theta_chief1, inc_chief, q1_chief, q2_chief, BigOmg_chief];
    AMap1 = ForwardMapping(OE_chief1, mu_Earth); % Linear mapping matrix
    
    rf = Xnomf(1,1:3)'; vf = Xnomf(1,4:6)';
    COEf = RVtoCOEs(rf,vf,mu_Earth);
    a2_chief = COEf(1); e2_chief = COEf(2); inc2_chief = COEf(3);
    BigOmg2_chief = COEf(4); LitOmg2_chief = COEf(5); f2 = COEf(6);
    q12_chief = e2_chief*cos(LitOmg2_chief); q22_chief = e2_chief*sin(LitOmg2_chief);
    theta2_chief = f2+LitOmg2_chief;
    OE_chief2 = [a2_chief, theta2_chief, inc2_chief, q12_chief, q22_chief, BigOmg2_chief];
    AMap2 = ForwardMapping(OE_chief2, mu_Earth); % Linear mapping matrix
    
    % Specify the final relative-orbital elements of the deputies in both orbits
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    dela = 0;
    dele = 25/(5*a_chief);
    deli = -50/(3*a_chief);
    delLitOmg = -2*pi*1E-6;
    delBigOmg = 0;
    delM = pi*1E-4;
    delq1_1 = dele*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_1 = [dela, deltheta, deli, delq1_1, delq2_1, delBigOmg]';
    
    
    dele = -20/(2*a_chief);
    deli = -25/(2*a_chief);
    delLitOmg = -5*pi*1E-4;
    delBigOmg = pi*1E-4;
    delM = 0;
    delq1_1 = dele*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_2 = [dela, deltheta, deli, delq1_1, delq2_1, delBigOmg]';
    
    
    dele = 20/(3*a_chief);
    deli = 50/(2*a_chief);
    delLitOmg = -2*pi*1E-6;
    delBigOmg = -5*pi*1E-6;
    delM = pi*1E-4;
    delq1_1 = dele*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_1 = dele*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    delCOE_3 = [dela, deltheta, deli, delq1_1, delq2_1, delBigOmg]';
    
    dele2 = 10/(8*a_chief);
    deli2 = 50/(4*a_chief);
    delBigOmg = -pi*1E-4;
    delM = 5*pi*1E-4;
    delLitOmg = pi*1E-3;
    deltheta = delLitOmg + delM; % Relative true latitude in rad
    
    delq1_2  = dele2*cos(BigOmg_chief) - e_chief*sin(BigOmg_chief)*delLitOmg;
    delq2_2  = dele2*sin(BigOmg_chief) + e_chief*cos(BigOmg_chief)*delLitOmg;
    delCOE_4 = [dela, deltheta, deli2, delq1_2, delq2_2, delBigOmg]';
    
    % Integrating the deputy initial and final trajectories
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    X01 = AMap1*delCOE_1;
    X02 = AMap1*delCOE_2;
    Xf1 = AMap2*delCOE_3;
    Xf2 = AMap2*delCOE_4;
    
    % q_t0 = load('q_t0.mat');  qr_chasser00 = q_t0.qr_chasser00;
    % q_tf = load('q_tf.mat');  qr_chasserf0 = q_tf.qr_chasserf0;
    qr_chasser00 = randi([1 9],[4 1]); %  [6; 4; 9; 1]; % [4; 3; 6; 8]; % 
    qr_chasserf0 = randi([1 9],[4 1]); %  [7; 5; 1; 8]; % [4; 2; 4; 3]; % 
    
    qr_chasser0 = qr_chasser00/norm(qr_chasser00);
    Omega_chaser0 = 5e-4*deg2rad(randi([1 9],[3 1])); % [4;8;4]  
    CN = Quaternion_to_DCM(qr_chasser0); Omega_chaser_B0 = CN*Omega_chaser0;
    X0 = [qr_chasser0; Omega_chaser0; X01];
    
    qr_chasserf = qr_chasserf0/norm(qr_chasserf0);
    Omega_chaserf = 3e-5*deg2rad(randi([1 9],[3 1])); % [4;4;2]
    CN = Quaternion_to_DCM(qr_chasserf); Omega_chaser_Bf = CN*Omega_chaserf;
    Xf = [qr_chasserf; Omega_chaserf; Xf1];
end
% -------------------------------------------------------------------------

%% ------- Nondimentionalized the inital condition of the problem ---------
for ll = 1
    options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-15);
    [~, Xrel0] = ode113(@(t,x)NonDimentionalized_NLode1(t,x,Target),TT1,X0(8:13),options);
    [~, Xrel1] = ode113(@(t,x)NonDimentionalized_NLode2(t,x,Target),TT2,Xf(8:13),options);
    [~, Xrel]  = ode113(@(t,x)DragSRPNonDimentionalized_FlatPlateOde(t,x,Target,Chasser),Tnorm,X0);
    for k = 1
        % TTime = linspace(0,SimTime,250);
        % for l = 1:250
        %     X(:,l) = mypdexic(TTime(l),X0,Xf,SimTime);
        % end
        % fh2 = figure;
        % hold on
        % plt1 = plot3(Xrel0(:,1),Xrel0(:,2),Xrel0(:,3),'Color',c5,'LineWidth',5);
        % plt2 = plot3(Xrel1(:,1),Xrel1(:,2),Xrel1(:,3),'Color',c3,'LineWidth',5);
        % plt3 = plot3(Xrel(:,8),Xrel(:,9),Xrel(:,10),'Color',c11,'LineWidth',5);
        % % plt7 = plot3(X(8,:),X(9,:),X(10,:),'Color',c9,'LineWidth',5);
        % plt4 = plot3(Xrel0(1,1),Xrel0(1,2),Xrel0(1,3),'o',...
        %     'LineWidth',2,...
        %     'MarkerEdgeColor',c3,...
        %     'MarkerFaceColor',c3,...
        %     'MarkerSize',10);
        % plt5 = plot3(Xrel1(1,1),Xrel1(1,2),Xrel1(1,3),'o',...
        %     'LineWidth',2,...
        %     'MarkerEdgeColor',c5,...
        %     'MarkerFaceColor',c5,...
        %     'MarkerSize',10);
        % plt6 = plot3(0,0,0,'o',...
        %     'LineWidth',2,...
        %     'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','k',...
        %     'MarkerSize',25);
        % grid on
        % view(-73,7)
        % grid on
        % grid minor
        % xlabel('X [km]')
        % ylabel('Y [km]')
        % zlabel('Z [km]')
        %
        % % hL = legend([plt1, plt2, plt3,plt4, plt5, plt6],...
        % % {'Initial Orbit','Final Orbit',...
        % % 'Uncontrolled Traj','Initial Condition',...
        % % 'Final Condition','Chief Spacecraft'},...
        % % 'AutoUpdate','off');
        % % saveLegendToImage(fh2, hL, './Figures/InitialConditionLegend')
        %
        % set(gca,'FontSize',40)
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
        % fh2.PaperPosition = [0 0 8580 7320]/res;
        % % print(fh2,'./Figures/3DInitialCondition','-dpng',sprintf('-r%d',res))
    end
end
% -------------------------------------------------------------------------

%% ------------------------ AutoDiff using CasADi -------------------------
for ll = 1
    addpath('../casadi-osx-matlabR2015a-v3.5.5')
    import casadi.*
    CasADiopt = struct('main', true, 'mex', true);
    t         = SX.sym('t');
    Xaug      = SX.sym('q1',13);  % Xaug  = [q1; omega1; rho1];
    dXaug     = SX.sym('q2',13);  % dXaug = [dq1;domega1;drho1];
    XChief    = SX.sym('XChief',6);
    Penalty1  = SX.sym('Penalty1');
    Penalty2  = SX.sym('Penalty2');
    Penalty3  = SX.sym('Penalty3');
    
    Ichaser_inv = 1./(Chasser.Moment_Of_Inertia_Calculator().');
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
            F = eye(N); F(5:7,5:7) = diag(Ichaser_inv); F_Aug = F;
            % penalty matrix (D matrix in the paper)
            D = diag([Penalty1*ones(1,4), ones(1,3), Penalty2*ones(1,3), Penalty3*ones(1,3)]);
            D_Aug = D;
        else
            F_Aug = blkdiag(F_Aug,F);
            D_Aug = blkdiag(D_Aug,D);
        end
    end
    % Defining the drift vector field
    drift = DragSRPNonDimentionalized_FlatPlateOdeCasADi(t,Xaug,Target,Chasser,XChief);
    
    % -------------------- Metric and curve Length ----------------------------
    % Defining the metric
    G = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);
    % actuated curve length
    L = Function('L',{Xaug,dXaug,Penalty1,Penalty2,Penalty3,t,XChief},...
        {(dXaug - drift).' * G * (dXaug - drift)},...
        {'X','dX','k1','k2','k3','t','XChief'},...
        {'CurveLength'});
    dLdx = L.jacobian();
    fdrift = Function('fdrift', {t,Xaug,XChief},...
        {drift});
    
    % Rewrite the functions in c files
    L.generate('L.c',CasADiopt);
    dLdx.generate('dLdx.c',CasADiopt);
    fdrift.generate('fdrift.c',CasADiopt);
    
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
%% --------------------------- AGHF parameters ----------------------------
for ll = 1
    clc
    smax = 5e6; % Need typo find an analytic way to solve for this value
    % # of grids in t
    tgrids = 250;
    % # of grids in s
    sgrids = 250;
    % # of grids of integration for final trajectory extraction
    intgrids = 1e5;
    % # of inputs
    N = length(X0);
    % penalty value (the lambda in paper)
    Penalty1 = 1e-3; Penalty2 = 5e1; Penalty3 = 5e2;
    % tolerance of the integrator
    % opts = odeset('RelTol',2.22045e-9,'AbsTol',2.22045e-14,'NormControl','on');
    opts = odeset('RelTol',2.22045e-9);
    % opts = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-14);
    % Setting the boundary condtions and the intial condition
    tmax = smax; xpoints = tgrids; tpoints = sgrids; m = 0;
    t = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of s interval in log scale % linspace(0,tmax,tpoints); % 
    T = SimTime; % motion duration
    x = linspace(0,T,xpoints); % [0 logspace(-4,log10(T),xpoints-1)];
    xnew = linspace(0,T,intgrids);
end
% -------------------------------------------------------------------------

%% ------------------------- Solving the AGHF -----------------------------
for ll = 1
    % solve for trajectory, see "AGHF" function file for implementation details
    clc
    tic;
    disp('Running PDEPE...');
    system('caffeinate -dims &');
    sol = pdepe(m,@(x,t,u,DuDx) mypdexpde(x,t,u,DuDx,Penalty1,Penalty2,Penalty3,N),...
                  @(x) mypdexic(x,X0,Xf,T),...
                  @(xl,ul,xr,ur,t) mypdexbc(xl,ul,xr,ur,t,X0,Xf,N),...
                  x,t,opts); % The solution "sol" is of form sol(t,x,i)
    system('killall caffeinate');
    toc;
    % save q_t0.mat qr_chasser00 -v7.3
    % save q_tf.mat qr_chasserf0 -v7.3
    % save sol.mat sol -v7.3
    % Sol = load('sol.mat');
    % sol = Sol.sol;
    
    % number = '2028171774';
    % carrier = 'AT&T';
    % message = 'Your MatLab code is done running!';
    % send_text_message(number,carrier,message)
end
% -------------------------------------------------------------------------

%% ------------------ Calculating the action integral ---------------------
for ll = 1
    tic;
    disp('Computing the action integral...');
    % calculate the atuated curve length for each grid of s
    X_temp  = zeros(N,xpoints);
    dX_temp = zeros(N,xpoints);
    Xval    = zeros(N,Tsize);
    cost    = zeros(tpoints,1);
    
    for j = 1:tpoints
        for kk = 1:N
            [X_temp(kk,:), dX_temp(kk,:)] = pdeval(m,x,sol(j,:,kk),x);
            if j == tpoints
                Xval(kk,:)                  = spline(x,X_temp(kk,:),Tnorm);
                AgentNom.DesiredState(:,kk) = spline(x,X_temp(kk,:),xnew)';
                AgentNom.FDesired(:,kk)     = spline(x,dX_temp(kk,:),xnew)';
            end
        end
        
        for i = 1:xpoints
            X_chief = full(ChiefState(x(i)));
            DuDx    = dX_temp(:,i);
            u       = X_temp(:,i);
            LL      = full(L(u,DuDx,Penalty1,Penalty2,Penalty3,x(i),X_chief));
            cost(j) = cost(j) + LL;
            
        end
    end
    AgentNom.FDesiredInterpolant     = griddedInterpolant(xnew,AgentNom.FDesired,'spline');
    AgentNom.DesiredStateInterpolant = griddedInterpolant(xnew,AgentNom.DesiredState,'spline');
    toc;
end
% -------------------------------------------------------------------------

%% ------------------ Integrate the system's trajectory -------------------
for ll = 1
    disp('Integrating the resulting trajectory...');
    options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-15); PMat=1e1*diag([1 2 3]); kMat = 3e1; % Controller gains
    tic
    [~,X_ode45] = ode15s(@(t,Xaug)DragSRP_FlatPlateControlledOde(t,Xaug,Target,Chasser,kMat,PMat,AgentNom),Tnorm,X0,options);
    toc
    disp('Done!!!!!');
end
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% ######################## Plotting Results ###############################
% -------------------------------------------------------------------------

%% ----------- Plotting the attitude quaternions of the deputy -----------
for ll = 1
    % close all
    fh2 = figure;
    LaBeL1 = {'$q_{0}$','$q_{1}$','$q_{2}$','$q_{3}$'};
    for i = 1:4
        subplot(2,2,i)
        plt1 = plot(tspan,Xval(i,:),'Color',c6,'LineWidth',2.5);
        hold on
        plt2 = plot(tspan,X_ode45(:,i).','Color',c11,'LineWidth',3);
        grid on
        % grid minor
        ylabel(LaBeL1{i})
        set(gca,'FontSize',25)
    end
    hL = legend([plt1, plt2],...
        {'AGHF','Lyapunov Ctrl'},'AutoUpdate','off','Location', 'Best');
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
    % print(fh2,'./Figures/Quaternion_Tracking_6DoF','-dpng',sprintf('-r%d',res))
end

%% --------- Plotting the action integral across each iteratons ----------
for ll = 1
    fh2 = figure;
    loglog(t(1:sgrids),cost,'Color',c5,'LineWidth',2.5);
    title('Action integral across iteration')
    xlabel('Homotopy iteration')
    ylabel('Length of the actuated curve')
    grid on;
    
    % set(findall(fh2,'Units','pixels'),'Units','normalized');
    % fh2.Units = 'pixels';
    % fh2.OuterPosition = [10 10 700 600];
    % res = 500;
    % set(fh2,'PaperPositionMode','manual')
    % fh2.PaperUnits = 'inches';
    % fh2.PaperPosition = [0 0 5580 4320]/res;
    % print(fh2,'./Figures/ActionIntegral','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% ---------------- Plotting the resulting 3D trajectories ----------------
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
        idx = 1+N*(jj-1):N*jj;
        
        h0 = plot3(Xval(idx(8),:),Xval(idx(9),:),Xval(idx(10),:),...
            '-','Color',c4,'LineWidth',2.5);
        % plot3(Xrel(:,8),Xrel(:,9),Xrel(:,10),'Color',c11,'LineWidth',2.5);
        % scatter3(X_temp(idx(8),:),X_temp(idx(9),:),X_temp(idx(10),:),...
        %     '*k','LineWidth',2.5);
        plt0(:,jj) = plot3(X_ode45(:,idx(8)),X_ode45(:,idx(9)),X_ode45(:,idx(10)),...
            '--','Color',cMat0(jj,:),'LineWidth',2.5);
        
        plt1(:,jj) = plot3(Xrel0(1,idx(1)),Xrel0(1,idx(2)),Xrel0(1,idx(3)),'*',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat1(jj,:),...
            'MarkerFaceColor',cMat1(jj,:)',...
            'MarkerSize',8);
        
        plt2(:,jj) = plot3(Xrel1(1,idx(1)),Xrel1(1,idx(2)),Xrel1(1,idx(3)),'o',...
            'LineWidth',3,...
            'MarkerEdgeColor',cMat2(jj,:),...
            'MarkerFaceColor',cMat2(jj,:),...
            'MarkerSize',8);
        
        h5 = plot3(0,0,0,'bo',...
            'LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','k',...
            'MarkerSize',15);
        plot3(Xrel0(:,idx(1)),Xrel0(:,idx(2)),Xrel0(:,idx(3)),'Color',cMat1(jj,:),'DisplayName','Initial Traj');
        plot3(Xrel1(:,idx(1)),Xrel1(:,idx(2)),Xrel1(:,idx(3)),'Color',cMat2(jj,:));
        
        grid on
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        title('{\color{black} 3D Relative Trajectory}','interpreter','tex')
        view(-73,23)
    end
    
    % hL = legend([h5, plt1(end,1), plt2(end,1), h0, plt0(end,1)],...
    %     {'Chief Location','Agent1 Init Condition',...
    %     'Agent1 End Condition', 'Converged Solution', 'Integrated Solution (ODE45)'},'AutoUpdate','off');
    % newPosition = [0.6 0.2 0.1 0.1];
    % newUnits = 'normalized';
    % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    % set(gca,'FontSize',25)
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
    % % save figure
    % print(fh2,'./Figures/Integrated_Trajectory_6DoF2','-dpng',sprintf('-r%d',res))
    
end
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% ######################## Required Functions #############################
% -------------------------------------------------------------------------
function [c,f,s] = mypdexpde(x,t,u,DuDx,k1,k2,k3,N) % Define PDE; right-hand-side of AGHF

XChief = full(ChiefState(x));
LL = full(L(u,DuDx,k1,k2,k3,x,XChief));
CasADiresult = full(dLdx(u,DuDx,k1,k2,k3,x,XChief,LL))';
CasADir = [CasADiresult(1:13), CasADiresult(14:26)]; % [dL/dx, dL/dxdot]

f = CasADir(:,2);  % dL/dxdot
s = -CasADir(:,1); % -dL/dx
c = ones(N,1);

end

function u0 = mypdexic(x,X0,Xf,T) % Initial condition for AGHF
% Sinusoidal initial condition
freq1 = pi/2 + 2*pi;
freq2 = pi + 2*pi;
u0 = X0*cos(freq1*x/T) ...
     + Xf*sin(freq1*x/T) ...
     - (X0+Xf)*sin(freq2*x/T);
u0(1:4) = u0(1:4)/norm(u0(1:4)); % imposing the quaternion constraint
end

function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t,X0,Xf,N) % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end
% -------------------------------------------------------------------------

function [Xdot] = NonDimentionalized_NLode1(t,X,Target)
global mu_Earth Tc 
format long
ti = t*Tc;

tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
X_chief = full(ChiefState0(t));
r_target = X_chief(1:3);  v_target = X_chief(4:6); rt_norm = (r_target.'*r_target).^(1/2);
rho_bar = X(1:3,:); drho_bar = X(4:6,:);

% Redimensionalizing the variables in the ODE
rho = rho_bar; rho_prime = drho_bar;
rc  = [rt_norm; 0; 0] + rho; rc_norm =  (rc.'*rc).^(1/2);

% Calculating the respective accelerations
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(r_target)*v_target; h_norm = norm(h_vec);  u_acc = U_eci_Target;
eh = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
u_acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                -2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(u_acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Rotation_Effect = - 2*tilde(Omega)*rho_prime ...
    - tilde(Omega)*(tilde(Omega)*rho) ...
    - tilde(Omegadot)*rho;

% Nondimetional relative ODE of the deputy
Xdot = [drho_bar; (del_ag + Rotation_Effect)];
Xdot = Xdot*Tc;

end

function [Xdot] = NonDimentionalized_NLode2(t,X,Target)
global mu_Earth Tc 
format long
ti = t*Tc;

tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
X_chief = full(ChiefStatef(t));
r_target = X_chief(1:3);  v_target = X_chief(4:6); rt_norm = (r_target.'*r_target).^(1/2);
rho_bar = X(1:3,:); drho_bar = X(4:6,:);

% Redimensionalizing the variables in the ODE
rho = rho_bar; rho_prime = drho_bar;
rc  = [rt_norm; 0; 0] + rho; rc_norm =  (rc.'*rc).^(1/2);

% Calculating the respective accelerations
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(r_target)*v_target; h_norm = norm(h_vec);  u_acc = U_eci_Target;
eh = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
u_acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                -2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(u_acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% relative gravitationall acceleration and coriolis effects
del_ag =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
Rotation_Effect = - 2*tilde(Omega)*rho_prime ...
    - tilde(Omega)*(tilde(Omega)*rho) ...
    - tilde(Omegadot)*rho;

% Nondimetional relative ODE of the deputy
Xdot = [drho_bar; (del_ag + Rotation_Effect)];
Xdot = Xdot*Tc;

end

function [Xdot] = DragSRPNonDimentionalized_FlatPlateOde(t,X,Target,Chasser)
format long
global mu_Earth Tc
ti = t*Tc;
tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];

% Collecting relevant quantities
X_chief = full(ChiefState(t));
r_target = X_chief(1:3); v_target = X_chief(4:6);
rt_norm = norm(r_target);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

q_Chasser = X(1:4)/norm(X(1:4));
Omega_Chasser_body = X(5:7);
rho = X(8:10); rho_prime = X(11:13);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
r_chasser = TN\rc; % ECI position of the chasser

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
FDeputy_body = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser);

%***************************************%
% Rotational state differential equation
I_Chasser        = diag(Chasser.Moment_Of_Inertia_Calculator());
OmegaMat_Chasser = [0, -Omega_Chasser_body(1), -Omega_Chasser_body(2), -Omega_Chasser_body(3);
                    Omega_Chasser_body(1),  0,  Omega_Chasser_body(3), -Omega_Chasser_body(2);
                    Omega_Chasser_body(2), -Omega_Chasser_body(3),  0,  Omega_Chasser_body(1);
                    Omega_Chasser_body(3), Omega_Chasser_body(2), -Omega_Chasser_body(1),  0];

qdot      = 1/2*OmegaMat_Chasser*q_Chasser;
Omega_dot = I_Chasser\(-tilde(Omega_Chasser_body)*(I_Chasser*Omega_Chasser_body));

%***************************************%
% Translational state differential equation
BN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = BN\FDeputy_body(1:3);

% Computing the angular velocity and angular acceleration of the target frame
h_vec   = tilde(r_target)*v_target; h_norm = norm(h_vec); u_acc = U_eci_Target;
eh      = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; Xrel_norm = norm(Xrel);
eh_dot  = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
acc_dot = beta_chief*(Vrel/Xrel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/Xrel_norm^5);

    
Omega    = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                - 2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% Computing the relative acceleration
Del_ag  =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelUsrp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP
Rotation_Effect = - 2*tilde(Omega)*rho_prime ...
                  - tilde(Omega)*(tilde(Omega)*rho) ...
                  - tilde(Omegadot)*rho;

% Differential Drag acceleration
Udrag_Target = Drag_CannonBall(r_target,v_target,Target);         
v_chasser = v_target + TN\(rho_prime + cross(Omega,rho));
UDrag_Chasser = Drag(r_chasser,v_chasser,q_Chasser,Chasser);
DelDrag = TN*(UDrag_Chasser-Udrag_Target);

% Integrating the relative trajectory
Rho_dot = [rho_prime; DelDrag + DelUsrp + Del_ag + Rotation_Effect];

% Collecting the rates of change
Xdot = [qdot; Omega_dot; Rho_dot]*Tc;

end
% -------------------------------------------------------------------------

function [Xdot] = DragSRPNonDimentionalized_FlatPlateOdeCasADi(t,X,Target,Chasser,X_chief)
format long
global mu_Earth Tc
ti = t*Tc;

tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];

% Collecting relevant quantities
r_target = X_chief(1:3); v_target = X_chief(4:6);
rt_norm = norm(r_target);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

q_Chasser = X(1:4)/norm(X(1:4));
Omega_Chasser_body = X(5:7);
rho = X(8:10);  rho_prime = X(11:13);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
r_chasser = TN\rc; % ECI position of the chasser

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
F_body = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser);

%***************************************%
% Rotational state differential equation
I_Chasser            = diag(Chasser.Moment_Of_Inertia_Calculator());
Omega_matrix_Chasser = [0, -Omega_Chasser_body(1), -Omega_Chasser_body(2), -Omega_Chasser_body(3);
                        Omega_Chasser_body(1),  0,  Omega_Chasser_body(3), -Omega_Chasser_body(2);
                        Omega_Chasser_body(2), -Omega_Chasser_body(3),  0,  Omega_Chasser_body(1);
                        Omega_Chasser_body(3), Omega_Chasser_body(2), -Omega_Chasser_body(1), 0];

qdot      = 1/2*Omega_matrix_Chasser*q_Chasser;
Omega_dot = I_Chasser\(-tilde(Omega_Chasser_body)*(I_Chasser*Omega_Chasser_body));

%***************************************%
% Translational state differential equation
BN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = BN\F_body(1:3);

% Computing the angular velocity and angular acceleration of the target frame
h_vec   = tilde(r_target)*v_target; h_norm = norm(h_vec); u_acc = U_eci_Target;
eh      = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; r_rel_norm = norm(Xrel);
eh_dot  = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
acc_dot = beta_chief*(Vrel/r_rel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/r_rel_norm^5);

    
Omega    = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                - 2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% Computing the relative acceleration
del_ag  =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelUsrp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP
Rotation_Effect = - 2*tilde(Omega)*rho_prime ...
                  - tilde(Omega)*(tilde(Omega)*rho) ...
                  - tilde(Omegadot)*rho;


% Integrating the relative trajectory
Rho_dot = [rho_prime; DelUsrp + del_ag + Rotation_Effect];

% Collecting the rates of change
Xdot = [qdot; Omega_dot; Rho_dot]*Tc;

end

