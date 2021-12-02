%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Hermann Kaptui Sipowa %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This work was presented at the 2022 IEEE/Aerospace Conference %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------

clear
close all
start_up
format long e
clc
global mu_Earth tspan T_normalized mu_dot ...
       Num_Agent a_chief e_chief Tc
%% ====================================================================== %
% ======================================================================= %
% ######################## Problem Inialization ######################### %
% ======================================================================= %
% ======================================================================= %

%% -------------------- Setting the Chief's parameters --------------------
for i = 1
    c1  = rgb('DarkGreen');
    c2  = rgb('DarkSlateGray');
    c3  = rgb('Aqua');
    c4  = rgb('LightSlateGray');
    c5  = rgb('DarkBlue');
    c6  = rgb('Red');
    c7  = rgb('Purple');
    c8  = rgb('Bisque');
    c9  = rgb('Orange');
    c10 = rgb('DarkGray');
    c11 = rgb('Teal');
    c12 = rgb('MediumVioletRed');
    c13 = rgb('SaddleBrown');
    % ----------------- Calculating the boundaries conditions ----------------
    mu_Earth    = 3.986004415E5; % Earth gravitational parameter
    Num_Agent   = 6; % Number of agents in the system
    Num_Agent1  = 1; % Number of agents in the first final-orbit
    Num_Agent2  = 1; % Number of agents in the second final-orbit
    M = 3;
    N = 6;
    % --------------------- Specifying the chief's orbit ----------------------
    a_chief      = 1.42e4;   % Semi-major axis in Km
    e_chief      = 0.5;      % Eccentricity
    inc_chief    = 50;       % Inclination in deg0
    BigOmg_chief = 30;       % RAAN in deg
    LitOmg_chief = 10;       % AOP in deg
    M_chief1     = 0;        % Initial mean anomaly
    
    
    % ------------------ Integrating the chief's trajectory ------------------
    SimTime = 0.75; % Simulation time (ratio of the chief's period)
    Period  = (2*pi)*sqrt(a_chief^3/mu_Earth);
    IntTime = SimTime*Period;
    tspan   = linspace(0,IntTime,1e4);
    options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-17);
    mu_dot  = sqrt(mu_Earth/a_chief^3); % Chief's mean motion
    Rc  = a_chief*eye(3); Tc = Period; T_normalized = tspan/Tc;
    inc_chief = deg2rad(inc_chief); BigOmg_chief = deg2rad(BigOmg_chief);
    LitOmg_chief = deg2rad(LitOmg_chief); M_chief1 = deg2rad(M_chief1);
    COE1 = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief1];
    [Position_target1,Velocity_target1] = COEstoRV(COE1,mu_Earth);
    X0_Chief1 = [Position_target1; Velocity_target1];
    [~, Xnom] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),tspan,X0_Chief1,options);
    
    Ttraj0     = linspace(0,Period,1e4);
    Ttrajf     = linspace(IntTime,Period+IntTime,1e4); X0_Chief2  = Xnom(end,:)';
    [~, Xnom0] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),Ttraj0,X0_Chief1,options);
    [~, Xnomf] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),Ttrajf,X0_Chief2,options);
    ChiefMeanMotion = sqrt(mu_Earth/a_chief^3);
    
    
    
end

%% ------- Create MX files to interpolate the chief's states (Xnom) -------
for i = 1
    % addpath('../casadi-osx-matlabR2015a-v3.5.5')
    % import casadi.*
    % CasADiopts = struct('main', true, 'mex', true);
    % TT0 = Ttraj0/Tc; TTf = Ttrajf/Tc; tt = MX.sym('t');
    %
    % %=====================================================================%
    % if exist('ChiefState','file') == 3
    %     delete ChiefState.c ChiefState.mexmaci64
    % end
    % Chief_x  = casadi.interpolant('Chief_x','bspline',{T_normalized},Xnom(:,1)');
    % Chief_y  = casadi.interpolant('Chief_y','bspline',{T_normalized},Xnom(:,2)');
    % Chief_z  = casadi.interpolant('Chief_z','bspline',{T_normalized},Xnom(:,3)');
    % Chief_dx = casadi.interpolant('Chief_dx','bspline',{T_normalized},Xnom(:,4)');
    % Chief_dy = casadi.interpolant('Chief_dy','bspline',{T_normalized},Xnom(:,5)');
    % Chief_dz = casadi.interpolant('Chief_dz','bspline',{T_normalized},Xnom(:,6)');
    % Output   = [Chief_x(tt);Chief_y(tt);Chief_z(tt);...
    %             Chief_dx(tt);Chief_dy(tt);Chief_dz(tt)];
    % ChiefState = Function('ChiefState',{tt},...
    %              {Output},{'time'},{'Output'});
    % ChiefState.generate('ChiefState.c',CasADiopts);
    % mex ChiefState.c -largeArrayDims
    %
    %
    % %=====================================================================%
    % if exist('ChiefState0','file') == 3
    %     delete ChiefState0.c ChiefState0.mexmaci64
    % end
    % Chief_x0  = casadi.interpolant('Chief_x0','bspline',{TT0},Xnom0(:,1)');
    % Chief_y0  = casadi.interpolant('Chief_y0','bspline',{TT0},Xnom0(:,2)');
    % Chief_z0  = casadi.interpolant('Chief_z0','bspline',{TT0},Xnom0(:,3)');
    % Chief_dx0 = casadi.interpolant('Chief_dx0','bspline',{TT0},Xnom0(:,4)');
    % Chief_dy0 = casadi.interpolant('Chief_dy0','bspline',{TT0},Xnom0(:,5)');
    % Chief_dz0 = casadi.interpolant('Chief_dz0','bspline',{TT0},Xnom0(:,6)');
    % Output0   = [Chief_x0(tt);Chief_y0(tt);Chief_z0(tt);...
    %              Chief_dx0(tt);Chief_dy0(tt);Chief_dz0(tt)];
    % ChiefState0 = Function('ChiefState0',{tt},...
    %               {Output0},{'time'},{'Output0'});
    % ChiefState0.generate('ChiefState0.c',CasADiopts);
    % mex ChiefState0.c -largeArrayDims
    %
    %
    % %=====================================================================%
    % if exist('ChiefStatef','file') == 3
    %     delete ChiefStatef.c ChiefStatef.mexmaci64
    % end
    % Chief_xf  = casadi.interpolant('Chief_xf','bspline',{TTf},Xnomf(:,1)');
    % Chief_yf  = casadi.interpolant('Chief_yf','bspline',{TTf},Xnomf(:,2)');
    % Chief_zf  = casadi.interpolant('Chief_zf','bspline',{TTf},Xnomf(:,3)');
    % Chief_dxf = casadi.interpolant('Chief_dxf','bspline',{TTf},Xnomf(:,4)');
    % Chief_dyf = casadi.interpolant('Chief_dyf','bspline',{TTf},Xnomf(:,5)');
    % Chief_dzf = casadi.interpolant('Chief_dzf','bspline',{TTf},Xnomf(:,6)');
    % Outputf   = [Chief_xf(tt);Chief_yf(tt);Chief_zf(tt);...
    %              Chief_dxf(tt);Chief_dyf(tt);Chief_dzf(tt)];
    % ChiefStatef = Function('ChiefStatef',{tt},...
    %               {Outputf},{'time'},{'Outputf'});
    % ChiefStatef.generate('ChiefStatef.c',CasADiopts);
    % mex ChiefStatef.c -largeArrayDims
    %
    % %=====================================================================%
    % clear Chief_x Chief_y Chief_z Chief_dx Chief_dy Chief_dz Output ChiefState ...
    %     Chief_x0 Chief_y0 Chief_z0 Chief_dx0 Chief_dy0 Chief_dz0 Output0 ChiefState0 ...
    %     Chief_xf Chief_yf Chief_zf Chief_dxf Chief_dyf Chief_dzf Outputf ChiefStatef
end

%% --------------- Setting the deputy initial conditions ------------------
for i = 1
    TT0 = Ttraj0/Tc; TTf = Ttrajf/Tc;
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
    X0 = nan(Num_Agent*N,1); Xf = nan(Num_Agent*N,1);
    % Specify the final relative-orbital elements of the deputies in both orbits
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    for l = 1
        for k = 1:Num_Agent
            idx = 1+N*(k-1):N*k;
            dela = 0;
        
            if k == 1
                aa = 2;
            else
                aa = k;
            end
        
            dele1     = -4*k/(3*a_chief);
            deli1     = -((-1)^k)*10/(aa*a_chief);
            delLitOmg = k*pi*1E-5;
            delBigOmg = ((-1)^k)*2*k*pi*1E-5;
            delM      = ((-1)^k)*k*pi/5*1E-4;
            delq1_1   = dele1*cos(LitOmg_chief) - e_chief*sin(LitOmg_chief)*delLitOmg;
            delq2_1   = dele1*sin(LitOmg_chief) + e_chief*cos(LitOmg_chief)*delLitOmg;
            deltheta  = delLitOmg + delM; % Relative true latitude in rad
            delCOE_0  = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
            X0(idx,1) = AMap1*delCOE_0;
        
        
            dele1     = k/(2*a_chief);
            deli1     = ((-1)^k)*10/(aa*a_chief);
            delLitOmg = (k/3)*pi*1E-4;
            delBigOmg = (k/15)*pi*1E-4;
            delM      = pi*1E-6;
            delq1_1   = dele1*cos(LitOmg_chief) - e_chief*sin(LitOmg_chief)*delLitOmg;
            delq2_1   = dele1*sin(LitOmg_chief) + e_chief*cos(LitOmg_chief)*delLitOmg;
            deltheta  = delLitOmg + delM; % Relative true latitude in rad
            delCOE_f  = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
            Xf(idx,1) = AMap2*delCOE_f;
        end
        % X02 = X0; shuffle1 = randperm(N,N);
        % for k = 1:length(shuffle1)
        %     j = shuffle1(k); idx  = 1+N*(k-1):N*k; idx2 = 1+N*(j-1):N*j;
        %     X0(idx,1) = X02(idx2,1);
        % end
        Xf2 = Xf; shuffle2 = [6 4 3 1 2 5];% randperm(N,N);  [6 3 1 4 2 5];
        for k = 1:length(shuffle2)
            j = shuffle2(k); idx  = 1+N*(k-1):N*k; idx2 = 1+N*(j-1):N*j;
            Xf(idx,1) = Xf2(idx2,1);
        end

    end % Nominal example used in the IEEE paper
    
    [~, Traj1] = ode113(@(t,x)NonDimentionalized_NLode0(t,x),TT0,X0,options);
    [~, Traj2] = ode113(@(t,x)NonDimentionalized_NLodef(t,x),TTf,Xf,options);
    AgentNom = struct('Xrel0',{Traj1(:,1:6),  Traj1(:,7:12),...
                               Traj1(:,13:18),Traj1(:,19:24),...
                               Traj1(:,25:30),Traj1(:,31:36)},...
                      'Xrelf',{Traj2(:,1:6),  Traj2(:,7:12),...
                               Traj2(:,13:18),Traj2(:,19:24),...
                               Traj2(:,25:30),Traj2(:,31:36)});
    %====================================================================
    for jj = 1
        % close all
        % LinWdth = 3;
        % MkrSize = 10;
        % fh2 = figure();
        % hold on
        % axis equal
        % view(47,7) % view(76,6)
        % grid on
        % grid minor
        % xlabel('X [km]')
        % ylabel('Y [km]')
        % zlabel('Z [km]')
        % set(gca,'FontSize',20)
        % set(findall(fh2,'Units','pixels'),'Units','normalized');
        % set(fh2, 'visible', 'on')
        % fh2.Units = 'pixels';
        % fh2.OuterPosition = [10 10 700 600];
        % cMat1 = [c9;c11;c3;c5;c6;c1];
        % for k = 1:Num_Agent
        %     idx = 1+N*(k-1):N*k;
        %     plot3(Traj1(:,idx(1)),Traj1(:,idx(2)),Traj1(:,idx(3)),'-','Color',cMat1(k,:),'LineWidth',5);
        %     plot3(Traj2(:,idx(1)),Traj2(:,idx(2)),Traj2(:,idx(3)),'-','Color',cMat1(k,:),'LineWidth',5);
        % end
        %
        %
        % plot3(0,0,0,'bo',...
        %     'LineWidth',1,...
        %     'MarkerEdgeColor','k',...
        %     'MarkerFaceColor','k',...
        %     'MarkerSize',20);
        %
        % plot3(X0(1),X0(2),X0(3),'s',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(X0(7),X0(8),X0(9),'s',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(X0(13),X0(14),X0(15),'s',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(X0(19),X0(20),X0(21),'s',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(X0(25),X0(26),X0(27),'s',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(X0(31),X0(32),X0(33),'s',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c6,...
        %     'MarkerFaceColor',c6',...
        %     'MarkerSize',MkrSize);
        %
        %
        % plot3(Xf(1),Xf(2),Xf(3),'p',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(Xf(7),Xf(8),Xf(9),'p',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(Xf(13),Xf(14),Xf(15),'p',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(Xf(19),Xf(20),Xf(21),'p',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(Xf(25),Xf(26),Xf(27),'p',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
        % plot3(Xf(31),Xf(32),Xf(33),'p',...
        %     'LineWidth',LinWdth,...
        %     'MarkerEdgeColor',c11,...
        %     'MarkerFaceColor',c11',...
        %     'MarkerSize',MkrSize);
    end
end

%% ====================================================================== %
% ======================================================================= %
% ############### Nominal Paths (No Replanning Guidance) ################ %
% ======================================================================= %
% ======================================================================= %

%% ----------- AutoDiff using CasADi For The Nominal Trajectory -----------
for j = 1
    % addpath('../casadi-osx-matlabR2015a-v3.5.5')
    % import casadi.*
    % CasADiopts = struct('main', true, 'mex', true);
    % XChief = SX.sym('XC',6);  tt    = SX.sym('t'); Penalty = SX.sym('Penalty');
    % rho1 = SX.sym('rho1',6);  drho1 = SX.sym('drho1',6);
    % rho2 = SX.sym('rho2',6);  drho2 = SX.sym('drho2',6);
    % rho3 = SX.sym('rho3',6);  drho3 = SX.sym('drho3',6);
    % rho4 = SX.sym('rho4',6);  drho4 = SX.sym('drho4',6);
    % rho5 = SX.sym('rho5',6);  drho5 = SX.sym('drho5',6);
    % rho6 = SX.sym('rho6',6);  drho6 = SX.sym('drho6',6);
    %
    % Xaug  = [rho1;rho2;rho3;rho4;rho5;rho6];
    % dXaug = [drho1;drho2;drho3;drho4;drho5;drho6];
    %
    % M = 3; N = 6;
    % if exist('L','file') == 3
    %     delete L.c L.mexmaci64 ...
    %            dLdx.c dLdx.mexmaci64 ...
    %            fdrift.c fdrift.mexmaci64
    % end
    %
    % for i = 1:Num_Agent
    %     if i == 1
    %         % [Fc F], note: the F_bar matrix in paper
    %         F = eye(N); F_Aug = F;
    %         % penalty matrix (D matrix in the paper)
    %         D = diag([Penalty*ones(1,N-M) ones(1,M)]); D_Aug = D;
    %     else
    %         F_Aug = blkdiag(F_Aug,F);
    %         D_Aug = blkdiag(D_Aug,D);
    %     end
    % end
    %
    % drift = NonDimentionalized_NLodeCasADi(tt,Xaug,XChief); % Drift vector field
    % H = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);
    % B = []; % Barrier Function (State Constraint)
    % G = (sum(B)+1)*H; % Riemannian Metric
    % L = Function('L',{Xaug,dXaug,Penalty,tt,XChief},...
    %     {(dXaug - drift).' * G * (dXaug - drift)},...
    %     {'X','dX','k','t','XChief'},...
    %     {'CurveLength'}); % Curve Length
    % dLdx = L.jacobian();
    % fdrift = Function('fdrift',{tt,Xaug,XChief},{drift});
    %
    % % Rewrite the functions as .Mex files
    % L.generate('L.c',CasADiopts);
    % dLdx.generate('dLdx.c',CasADiopts);
    % fdrift.generate('fdrift.c',CasADiopts);
    %
    % % Generating Mex files from the c files
    % mex L.c -largeArrayDims
    % mex dLdx.c -largeArrayDims
    % mex fdrift.c -largeArrayDims
    %
    % clear L dLdx fdrift rho1 drho1 rho2 drho2
end

%% --------------------------- AGHF parameters ----------------------------
for ll = 1
    clc
    smax = 1e7; % Need tpo find an analytic way to solve for this value
    % # of grids in t
    tgrids = 250;
    % # of grids in s
    sgrids = 250;
    % # of grids of integration for final trajectory extraction
    intgrids = 1e5;
    % # of inputs
    M = 3;
    % penalty value (the lambda in paper)
    Penalty = 5e5;
    % tolerance of the integrator
    % opts = odeset('AbsTol',1e-14);
    opts = odeset('RelTol',2.22045e-8,'AbsTol',2.22045e-20);
    
    % Setting the boundary condtions and the intial condition
    tmax = smax; xpoints = tgrids; xpoints_new = intgrids; tpoints = sgrids;
    m = 0; t = [0 logspace(-2,log10(tmax),tpoints-1)]; % discretization of s interval in log scale
    T    = SimTime; % motion duration
    x    = linspace(0,T,xpoints);
    xnew = linspace(0,T,intgrids);
    Time = Tc*xnew/3600;
end

%% ------------------------- Solving the AGHF -----------------------------
for ll = 1
% solve for trajectory, see "AGHF" function file for implementation details
close all
clc
tic;
disp('solving PDE...');

% Solve AGHF
sol = pdepe(m,@(x,t,u,DuDx) mypdexpde(x,t,u,DuDx,Penalty,length(X0)),...
    @(x) mypdexic(x,X0,Xf,T),...
    @(xl,ul,xr,ur,t) mypdexbc(xl,ul,xr,ur,t,X0,Xf,length(X0)),...
    x,t); % The solution "sol" is of form sol(t,x,i)
toc;
end

%% -- Tracking the Computed Optimal Path (Lyapunov Tracking Controller) --
for ll = 1
    for k = 1:Num_Agent
        TimeExtended = [xnew TTf(2:end)];
        idx = 1+N*(k-1):N*k;
        for i = 1:6
            [X_End, dX_End] = pdeval(m,x,sol(end,:,idx(i)),x);
            AgentNom(k).FDesired(:,i)     = spline(x,dX_End,xnew)';
            AgentNom(k).DesiredState(:,i) = spline(x,X_End,xnew)';
        end
        StatesExtended = [AgentNom(k).DesiredState; AgentNom(k).Xrelf(2:end,:)];
        AgentNom(k).FDesiredInterpolant       = griddedInterpolant(xnew,AgentNom(k).FDesired,'spline');
        AgentNom(k).DesiredStateInterpolant   = griddedInterpolant(xnew,AgentNom(k).DesiredState,'spline');
        AgentNom(k).StatesExtendedInterpolant = griddedInterpolant(TimeExtended,StatesExtended,'spline');
        
        %=====================================================================%
        FMat = eye(N); DMat = diag([Penalty*ones(1,N-M) ones(1,M)]);
        GMat = (FMat.')^(-1)*DMat*FMat^(-1);
        % Controller's gains
        KMat = eye(3)*3E-3;
        PMat = eye(3)*1E-1;
        AgentID = k;
        [~, AgentNom(k).X_LyapTrack] = ...
            ode113(@(t,X)LyapunovControlledTrajectory(t,X,KMat,PMat,GMat,AgentID,AgentNom),xnew,X0(idx),options);
        AgentNom(k).DesiredState = AgentNom(k).X_LyapTrack;
        for j = 1:length(xnew)
            [AgentNom(k).Crtl(j,:), AgentNom(k).LuapCost(j,1), AgentNom(k).AGHFCost(j,1)] = ...
                LyapunovTrackingCtrl(xnew(j),AgentNom(k).X_LyapTrack(j,:).',KMat,PMat,GMat,AgentID,AgentNom);
        end
    end
end
for ll = 1
    % close all
    % cMat0 = [c6;c1;c3;c4;c7;c12];
    % cMat1 = [c5;c5;c5;c5;c5;c5];
    % cMat2 = [c11;c11;c11;c11;c11;c11];
    % fh2 = figure;
    % view(68,10)
    % hold on
    % plt00 = zeros(Num_Agent);
    % plt0 = zeros(Num_Agent);
    % plt1 = zeros(Num_Agent);
    % plt2 = zeros(Num_Agent);
    %
    % for k = 1:Num_Agent
    %     AgentID = 1+N*(k-1):N*k;
    %     plt0(:,k) = plot3(AgentNom(k).DesiredState(:,1),AgentNom(k).DesiredState(:,2),AgentNom(k).DesiredState(:,3),...
    %         '-.','Color',cMat0(k,:),'LineWidth',5);
    %
    %     plt1(:,k) = plot3(X0(AgentID(1)),X0(AgentID(2)),X0(AgentID(3)),'*',...
    %         'LineWidth',5,...
    %         'MarkerEdgeColor',cMat1(k,:),...
    %         'MarkerFaceColor',cMat1(k,:)',...
    %         'MarkerSize',12);
    %
    %     plt2(:,k) = plot3(Xf(AgentID(1)),Xf(AgentID(2)),Xf(AgentID(3)),'o',...
    %         'LineWidth',5,...
    %         'MarkerEdgeColor',cMat2(k,:),...
    %         'MarkerFaceColor',cMat2(k,:),...
    %         'MarkerSize',12);
    %
    %     h5 = plot3(0,0,0,'bo',...
    %         'LineWidth',5,...
    %         'MarkerEdgeColor','k',...
    %         'MarkerFaceColor','k',...
    %         'MarkerSize',20);
    % grid on
    % grid minor
    % xlabel('X [km]')
    % ylabel('Y [km]')
    % zlabel('Z [km]')
    %
    % set(gca,'FontSize',20)
    % % set all units inside figure to normalized so that everything is scaling accordingly
    % set(findall(fh2,'Units','pixels'),'Units','normalized');
    % % do show figure on screen
    % set(fh2, 'visible', 'on')
    % % set figure units to pixels & adjust figure size
    % fh2.Units = 'pixels';
    % fh2.OuterPosition = [10 10 700 600];%[10 10 1400 1000];
    % % define resolution figure to be saved in dpi
    % res = 500;
    % % recalculate figure size to be saved
    % set(fh2,'PaperPositionMode','manual')
    % fh2.PaperUnits = 'inches';
    % fh2.PaperPosition = [0 0 8580 6320]/res;
    % end
end

%% ====================================================================== %
% ======================================================================= %
% ###################### Replanning Guidance ############################ %
% ======================================================================= %
% ======================================================================= %

%% ---------- AutoDiff using CasADi For The Replanning Algorithm ----------
for j = 1
    % addpath('../casadi-osx-matlabR2015a-v3.5.5')
    % import casadi.*
    % CasADiopts = struct('main', true, 'mex', true);
    % XChief = SX.sym('XC',6); tt = SX.sym('tt'); Penalty = SX.sym('Penalty');
    % Lr1 = SX.zeros(1); Lr2 = SX.zeros(1); Lr3 = SX.zeros(1);
    % Lr4 = SX.zeros(1); Lr5 = SX.zeros(1); Lr6 = SX.zeros(1);
    % Lr7 = SX.zeros(1); Lr8 = SX.zeros(1); Lr9 = SX.zeros(1);
    %
    % if exist('Lreplan','file') == 3
    %     delete Lreplan.c Lreplan.mexmaci64 ...
    %            Lreplandx.c Lreplandx.mexmaci64
    % end
    % if exist('fdriftreplan','file') == 3
    %     delete fdriftreplan.c fdriftreplan.mexmaci64
    % end
    % rho   = SX.sym('rho',6);   drho  = SX.sym('drho',6);
    % Obst1 = SX.sym('Obst1',3); Obst2 = SX.sym('Obst2',3);
    % Obst3 = SX.sym('Obst3',3); Obst4 = SX.sym('Obst4',3);
    % Obst5 = SX.sym('Obst5',3); Obst6 = SX.sym('Obst6',3);
    % Obst7 = SX.sym('Obst7',3); Obst8 = SX.sym('Obst8',3);
    % Obst9 = SX.sym('Obst9',3);
    % ObstaclesStates = [Obst1 Obst2 Obst3 Obst4 Obst5 Obst6 Obst7 Obst8 Obst9];
    %
    % % Defining the metric, before adding the state constraint barrier functions
    % F = eye(N); D = diag([Penalty*ones(1,N-M) ones(1,M)]); H = (F.')^(-1)*D*F^(-1);
    % % Defining the drift vector field
    % drift = NonDimentionalized_NLodeCasADi(tt,rho,XChief);
    %
    % % ------------------- state constraint barrier function ------------------
    % B  = SX.zeros(1); Num_Obstacles = Num_Agent-1;
    % kb = SX.sym('Kb',Num_Agent,2); pb = SX.sym('Pb',Num_Agent,2);
    %
    % for i = 1:Num_Obstacles
    %     delR = rho(1:3) - ObstaclesStates(:,i);
    %     B = B + kb(i,1)*exp( -(1/2)*(norm(delR)/pb(i,1))^2 );
    %     G = (B+1)*H;
    %
    %     if i==1
    %         Lr1 = (drho - drift).' * G * (drho - drift);
    %     elseif i==2
    %         Lr2 = (drho - drift).' * G * (drho - drift);
    %     elseif i==3
    %         Lr3 = (drho - drift).' * G * (drho - drift);
    %     elseif i==4
    %         Lr4 = (drho - drift).' * G * (drho - drift);
    %     elseif i==5
    %         Lr5 = (drho - drift).' * G * (drho - drift);
    %     elseif i==6
    %         Lr6 = (drho - drift).' * G * (drho - drift);
    %     elseif i==7
    %         Lr7 = (drho - drift).' * G * (drho - drift);
    %     elseif i==8
    %         Lr8 = (drho - drift).' * G * (drho - drift);
    %     elseif i==9
    %         Lr9 = (drho - drift).' * G * (drho - drift);
    %     end
    % end
    %
    % Lreplan = Function('Lreplan',{rho,drho,Penalty,tt,XChief,ObstaclesStates,kb,pb},...
    %     {Lr1,Lr2,Lr3,Lr4,Lr5,Lr6,Lr7,Lr8,Lr9},...
    %     {'X','dX','k','t','XChief','ObstaclesStates','Kb','Pb'},...
    %     {'Lr1','Lr2','Lr3','Lr4','Lr5','Lr6','Lr7','Lr8','Lr9'});
    % Lreplandx = Lreplan.jacobian();
    % Lreplan.generate('Lreplan.c',CasADiopts);
    % Lreplandx.generate('Lreplandx.c',CasADiopts);
    % mex Lreplan.c -largeArrayDims
    % mex Lreplandx.c -largeArrayDims
    %
    % fdriftreplan = Function('fdriftreplan',{tt,rho,XChief},{drift});
    % fdriftreplan.generate('fdriftreplan.c',CasADiopts);
    % mex fdriftreplan.c -largeArrayDims
    % % clear fd Lreplan Lreplandx
end

%% -------------------- Receding Horizon Path Planner ---------------------
for jk = 1
    clc
    tmax  = 1e7;
    t     = [0 logspace(-2,log10(tmax),tpoints-1)];
    FMat2 = eye(N); 
    DMat2 = diag([Penalty*ones(1,N-M) ones(1,M)]);
    HMat2 = (FMat2.')^(-1)*DMat2*FMat2^(-1);
    mins  = 60/Tc;     % Mapping from minutes to nondimesional time
    freq  = 10*mins;   % Guidance's update frquency
    Tbeta = 5*mins;    % Time to go to reach the edge of replanning set
    Tpsi  = 30*mins;   % Time to go to reach the edge of the Keep-out set
    ConeAngle  = 30*pi/180; % Conning angle used to define the keep-out zone
    Priority   = 1:Num_Agent;     % Array of the priority numner of each agent
    Separation = nan(Num_Agent); % Matrix of the separation distance between agents
    Partition  = freq*(0:floor(xnew(end)/freq)).';
    Penalty    = 5e3;
    kbMat      = 1e5*Priority;
    pb         = nan(Num_Agent,2);
    kb         = nan(Num_Agent,2);
    RelativeAngle = nan(Num_Agent);
    Occurence_idx = nan(length(Partition),1);
    i1 = 1; i2 = 1; i3 = 1; i4 = 1; i5 = 1;
    for i = 1:Num_Agent
        AgentNom(i).DesiredState = AgentNom(i).X_LyapTrack;
        AgentNom(i).DesiredStateInterpolant = griddedInterpolant(xnew,AgentNom(i).DesiredState,'spline');
        AgentNom(i).CrtlNew = AgentNom(i).Crtl;
    end
    for i = 1:length(Partition)
        Occurence_idx(i) = find(xnew>=Partition(i), 1, 'first');
    end
    for i = 1:length(Partition)
        idx = Occurence_idx(i);
        tk  = xnew(idx);
        KeepOutRadius = ReachableSetRadius(tk,Tpsi,AgentNom,Num_Agent); % Radius of the reachable sets
        
        % Computing the separation between agents in the system
        for j = 1:Num_Agent
            for k = 1:Num_Agent
                if j~=k
                    vel  = AgentNom(j).DesiredState(idx,4:6);
                    delR = AgentNom(k).DesiredState(idx,1:3)...
                        - AgentNom(j).DesiredState(idx,1:3);
                    Separation(j,k)    = norm(delR);
                    RelativeAngle(j,k) = acos(dot(delR,vel)/(norm(delR)*norm(vel)));
                end
            end
        end
        
        for j = 1:Num_Agent
            % Check if agent i needs to replan its path
            OtherAgentidx = 1:Num_Agent~=j;
            ObstaclesList = find(Separation(j,:) <= KeepOutRadius(j) & Priority(j) < Priority...
                & RelativeAngle(j,:)<ConeAngle);
            if ~isempty(ObstaclesList)
                if j == 1
                    AgentNom(j).IterationIndex(i1) = idx;
                    AgentNom(j).KeepOutRadius(i1)  = KeepOutRadius(j);
                    i1 = i1+1;
                elseif j == 2
                    AgentNom(j).IterationIndex(i2) = idx;
                    AgentNom(j).KeepOutRadius(i2)  = KeepOutRadius(j);
                    i2 = i2+1;
                elseif j == 3
                    AgentNom(j).IterationIndex(i3) = idx;
                    AgentNom(j).KeepOutRadius(i3)  = KeepOutRadius(j);
                    i3 = i3+1;
                elseif j == 4
                    AgentNom(j).IterationIndex(i4) = idx;
                    AgentNom(j).KeepOutRadius(i4)  = KeepOutRadius(j);
                    i4 = i4+1;
                elseif j == 5
                    AgentNom(j).IterationIndex(i5) = idx;
                    AgentNom(j).KeepOutRadius(i5)  = KeepOutRadius(j);
                    i5 = i5+1;
                end
                
                
                % Identifying the moving ostacles
                treplan = tk;%+Tbeta;
                xreplan = linspace(treplan,T,50);
                
                %=================== Predict the states of the obstacle ===================%
                for ijk = 1:length(ObstaclesList)
                    X0_replan = AgentNom(ObstaclesList(ijk)).DesiredState(idx,:)';
                    [~, Xrel] = ode113(@(t,x)NonDimentionalized_NLode(t,x),xreplan,X0_replan,options);
                    % AgentNom(ObstaclesList(ijk)).ObstaclePredictionInterpolant = griddedInterpolant(xreplan,Xrel,'spline');
                    AgentNom(ObstaclesList(ijk)).ObstaclePredictionInterpolant = ...
                        griddedInterpolant(xnew(idx:end),AgentNom(ijk).DesiredState(idx:end,:),'spline');
                end
                
                %=========================== Replanning sequence ===========================%
                AgentX0 = AgentNom(j).DesiredState(idx,:)'; AgentXf = Xf(1+N*(j-1):N*j);
                pb(1:length(ObstaclesList),:) = [Separation(ObstaclesList,j)/2 ObstaclesList']; 
                kb(1:length(ObstaclesList),:) = [kbMat(ObstaclesList)' ObstaclesList'];
                tic
                NewPath = pdepe(m,@(x,t,u,DuDx) ReplanningPDE(x,t,u,DuDx,Penalty,N,ObstaclesList,kb,pb,AgentNom),...
                    @(x) ReplanningIC(x,AgentNom,j,tk,AgentX0,AgentXf,T),...
                    @(xl,ul,xr,ur,t) ReplanningBCs(xl,ul,xr,ur,t,AgentX0,AgentXf,N),...
                    xreplan,t);
                toc
                
                %=================== Updating the trajectory of the deputies ===================%
                idxPlanned = idx;
                for ii = 1:6
                    [X_End, dX_End] = pdeval(m,xreplan,NewPath(end,:,ii),xreplan);
                    AgentNom(j).FDesired(idx:end,ii)     = spline(xreplan,dX_End,xnew(idx:end))';
                    AgentNom(j).DesiredState(idx:end,ii) = spline(xreplan,X_End,xnew(idx:end))';
                end
                AgentNom(j).FDesiredInterpolant     = griddedInterpolant(xnew,AgentNom(j).FDesired,'spline');
                AgentNom(j).DesiredStateInterpolant = griddedInterpolant(xnew,AgentNom(j).DesiredState,'spline');
                
                %=================== Updating the required control for path tracking ===================%
                for lk = idx:length(xnew)
                    rhoAgent = AgentNom(j).DesiredStateInterpolant(xnew(lk))'; B = [];
                    for kij = 1:length(ObstaclesList)
                        l = ObstaclesList(kij);
                        rhoObstacle = AgentNom(l).DesiredStateInterpolant(xnew(lk))';
                        delR = rhoAgent(1:3) - rhoObstacle(1:3);
                        B = kb(1,1)*exp( -(1/2)*(norm(delR)/pb(1,1))^2 );
                    end
                    GMat2 = (B+1)*HMat2;
                    [AgentNom(j).CrtlNew(lk,:),~,AgentNom(j).AGHFCostNew(lk,1)] = ...
                        LyapunovTrackingCtrl(xnew(lk),AgentNom(j).DesiredState(lk,:).',KMat,PMat,GMat2,j,AgentNom);
                end
            end
        end
    end

    for k = 1:Num_Agent
        NumbArrows = 10;
        tk = 0;
        M_max = M_chief1 + T*Tc*ChiefMeanMotion;
        dt = M_max/(NumbArrows*Tc*ChiefMeanMotion);
        ijk = 1;
        for ii = 1:NumbArrows-4
            tk2 = tk+0.35/(1.5*NumbArrows - ii);
            r1_base = AgentNom(k).DesiredStateInterpolant(tk);
            r1_head = AgentNom(k).DesiredStateInterpolant(tk2);
            Arrow = (r1_head-r1_base)/norm(r1_head-r1_base);
            if ijk == ii-2
                AgentNom(k).ArrowsBase(ijk,:) = r1_base(1:3);
                AgentNom(k).ArrowsHead(ijk,:) = r1_base(1:3) + Arrow(1:3);
                ijk = ijk+1;
            end
            tk = tk+dt;
        end
        AgentNom(k).TotalCtrl_Int = 1e6*simpsons(vecnorm(AgentNom(k).Crtl,2,2),0,0.75,[]);
        AgentNom(k).TotalCtrl_New = 1e6*simpsons(vecnorm(AgentNom(k).CrtlNew,2,2),0,0.75,[]);
        AgentNom(k).CrtPerCentChange = 100*(AgentNom(k).TotalCtrl_New - AgentNom(k).TotalCtrl_Int)/AgentNom(k).TotalCtrl_Int;
    end
end

%% ------------------------------------------------------------------------
% ######################## Plotting Results ###############################
% -------------------------------------------------------------------------
for jk = 1
    % clc
    close all
    cMat0 = [c11;c11;c3];
    cMat1 = [c9;c5;c3];
    cMat2 = [c3;c1];
    fh2 = figure;
    for jj = 1:Num_Agent
        % fh2 = figure;
        linesize = 5;
        % axis equal
        view(71,10)
        grid on
        % grid minor
        xlabel('X [km]')
        ylabel('Y [km]')
        zlabel('Z [km]')
        % title(sprintf('Agent %d',jj))
        hold on
        
        plt3 = plot3(AgentNom(jj).DesiredState(:,1),...
            AgentNom(jj).DesiredState(:,2),AgentNom(jj).DesiredState(:,3),...
            '-.','Color',c9,'LineWidth',linesize);
        
        plt2 = plot3(AgentNom(jj).X_LyapTrack(:,1),...
            AgentNom(jj).X_LyapTrack(:,2),AgentNom(jj).X_LyapTrack(:,3),...
            '-','Color',c11,'LineWidth',linesize);
        
        plot3(AgentNom(jj).X_LyapTrack(1,1),...
            AgentNom(jj).X_LyapTrack(1,2),AgentNom(jj).X_LyapTrack(1,3),...
            's','LineWidth',4,'MarkerEdgeColor',cMat0(2,:),...
            'MarkerFaceColor',c11,'MarkerSize',15);
        
        plot3(AgentNom(jj).X_LyapTrack(end,1),...
            AgentNom(jj).X_LyapTrack(end,2),AgentNom(jj).X_LyapTrack(end,3),...
            'p','LineWidth',4,'MarkerEdgeColor',cMat0(2,:),...
            'MarkerFaceColor',c11,'MarkerSize',20);
        
        
        
        % vel = AgentNom(jj).DesiredState(idx,4:6)'; b3 = vel/norm(vel); zhat = [0;0;1];
        % b2 = cross(b3,zhat); b2 = b2/norm(b2);
        % b1 = cross(b2,b3);   b1 = b1/norm(b1);
        % R = [b1 b2 b3];
        % MotionCone(KeepOutRadius(1),ConeAngle,R,AgentNom(1).DesiredState(idx,1:3),c11)
        %
        % p1 = AgentNom(1).DesiredState(idx,1:3)';
        % p2 = AgentNom(1).DesiredState(idx,1:3)' + 1500*AgentNom(1).DesiredState(idx,4:6)';
        % arrow3d(p1',p2',[c11;c11],[0.5,0.3])
        
        plt8 = arrow3d(AgentNom(jj).ArrowsBase,AgentNom(jj).ArrowsHead,[c6;c6],[0.45,0.5]);
        % plt9 = arrow3d(Agent2ArrowsBase,Agent2ArrowsHead,[c6;c6],[0.5,0.5]);
        % arrow3d(AgentNom(2).X_LyapTrack(idx,1:3),AgentNom(1).X_LyapTrack(idx,1:3),[c12;c12],[0.15,0.25]);
        % r = KeepOutRadius(2);
        % [x_sphere, y_sphere, z_sphere] = sphere(100);
        % h = surf(CheckLocation(1)+r*x_sphere, CheckLocation(2)+r*y_sphere, CheckLocation(3)+r*z_sphere);
        % set(h,'FaceAlpha', 0.09, 'EdgeColor', 'b', 'EdgeAlpha',.05, 'FaceColor', 'b');
        
        % hL = legend([plt2, plt3],...
        %     {'Initial Paths',...
        %     'Replanned Paths'},...
        %     'AutoUpdate','off');
        % hL.FontSize = 40;
        % newPosition = [0.25 0.25 0.1 0.1];
        % newUnits = 'normalized';
        % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
        % saveLegendToImage(fh2, hL, './Figures/3DLegend')
        
        set(gca,'FontSize',40)
        % set all units inside figure to normalized so that everything is scaling accordingly
        set(findall(fh2,'Units','pixels'),'Units','normalized');
        % do show figure on screen
        set(fh2, 'visible', 'on')
        % set figure units to pixels & adjust figure size
        fh2.Units = 'pixels';
        fh2.OuterPosition = [10 10 1000 800]; % [10 10 1400 1000];
    end
    % % define resolution figure to be saved in dpi
    % res = 500;
    % % recalculate figure size to be saved
    % set(fh2,'PaperPositionMode','manual')
    % fh2.PaperUnits = 'inches';
    % fh2.PaperPosition = [0 0 8580 6320]/res;
    % print(fh2,'./Figures/3DPlotSixAgents','-dpng',sprintf('-r%d',res))
    
    %======================================================================%
    % close all
    fh = figure;
    X = categorical({'Agent1','Agent2','Agent3','Agent4','Agent5','Agent6'});
    X = reordercats(X,{'Agent1','Agent2','Agent3','Agent4','Agent5','Agent6'});
    Y = nan(Num_Agent,2); ii = 1;
    for k = 1:Num_Agent
        Y(k,:) = [AgentNom(k).TotalCtrl_Int AgentNom(k).TotalCtrl_New];
        if isempty(AgentNom(k).AGHFCostNew)
            noreplan(ii) = k;
            ii = ii+1;
        end
    end
    noreplan(3) = 1;
    hb = bar(X,Y,1);
    hb(1).FaceColor = '#826E91'; hb(1).EdgeColor = 'none';
    hb(2).FaceColor = '#91826E';  hb(2).EdgeColor = 'none';
    hold on
     
    xtips2 = hb(2).XEndPoints;
    ytips2 = hb(2).YEndPoints;
    labels2 = string(hb(2).YData); labels2(noreplan) = '';
    text(xtips2,ytips2+1e-2,labels2,'HorizontalAlignment','left',...
        'VerticalAlignment','middle','FontSize',30,'Rotation',90)
    xtips1 = hb(1).XEndPoints; xtips1(noreplan) = xtips1(noreplan)+1.6e-1; 
    ytips1 = hb(1).YEndPoints; ytips1(noreplan) = ytips1(noreplan)+2e-2;
    labels1 = string(hb(1).YData); % labels1(end) = 'Initial Path';
    text(xtips1,ytips1+5e-3,labels1,'HorizontalAlignment','left',...
        'VerticalAlignment','middle','FontSize',30,'Rotation',90)
    
    for iji = 1:length(noreplan)
        iidx = noreplan(iji);
        bX = [xtips1(iidx)-2.0e-1, xtips1(iidx)-2.0e-1, xtips1(iidx)-1.2e-1,...
            xtips1(iidx)+0.05e-1, xtips2(iidx)-0.2e-1, ...
            xtips2(iidx)+.5e-1, xtips2(iidx)+.5e-1];
        bY = [ytips2(iidx), ytips2(iidx)+1e-2, ytips2(iidx)+1e-2...
            ytips2(iidx)+2e-2, ytips2(iidx)+1e-2...
            ytips2(iidx)+1e-2, ytips2(iidx)];
        plot(bX,bY, 'LineWidth', 3, 'Color', 'k');
    end
    
    xtips2(noreplan) = xtips2(noreplan)-1.5e-1;
    ytips2(noreplan) = ytips2(noreplan)-ytips2(noreplan)./(1.1);
    labels2(noreplan) = 'No Replanning';
    text(xtips2(noreplan),ytips2(noreplan),labels2(noreplan),'HorizontalAlignment','left',...
        'VerticalAlignment','middle','FontSize',40,'Rotation',90,'Color',c3,...
        'FontSize', 50, 'FontWeight', 'bold')
    
    ylim([0 max([ytips1 ytips2]+1e-1)])
    set(gca,'ytick',[])
    ylabel('$||\vec{u}||$[mm/s$^2$]')
    
    hL = legend([hb(1), hb(2)],...
        {'Initial Paths$~~~~~~~~~~~~~~~$',...
        'Replanned Paths'},...
        'AutoUpdate','off');
    hL.FontSize = 40;
    newPosition = [0 0 0.5 0.5];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',2,'Location','southoutside');
    
    
    
    set(gca,'FontSize',40)
    % set all units inside figure to normalized so that everything is scaling accordingly
    set(findall(fh,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh.Units = 'pixels';
    fh.OuterPosition = [10 10 600 500]; % [10 10 1400 1000];
    
    % res = 500;
    % set(fh,'PaperPositionMode','manual')
    % fh.PaperUnits = 'inches';
    % fh.PaperPosition = [0 0 8580 6320]/res;
    % print(fh,'./Figures/CtrlComparison','-dpng',sprintf('-r%d',res))
    
    
    % %======================================================================%
    % % close all
    % cMat0 = [c11;c13;c5];
    % cMat1 = [c9;c5;c5];
    % cMat2 = [c3;c1];
    % fh2 = figure;
    % textcase = [4,5,6];
    % plt1 = zeros(3);
    % for ij = 1:3%1:Num_Agent
    %     % fh2 = figure;
    %     jj = textcase(ij);
    %     linesize = 5;
    %     % axis equal
    %     view(100,10)
    %     hold on
    %     grid on
    %     grid minor
    %     xlabel('X [km]')
    %     ylabel('Y [km]')
    %     zlabel('Z [km]')
    %     % title(sprintf('Agent %d',jj))
    %     if ij==1
    %     plt0 = plot3(AgentNom(jj).DesiredState(:,1),...
    %         AgentNom(jj).DesiredState(:,2),AgentNom(jj).DesiredState(:,3),...
    %         '-.','Color',cMat1(jj-3,:),'LineWidth',linesize);
    %     end
    %     plt1(:,ij) = plot3(AgentNom(jj).X_LyapTrack(:,1),...
    %         AgentNom(jj).X_LyapTrack(:,2),AgentNom(jj).X_LyapTrack(:,3),...
    %         '-','Color',cMat0(ij,:),'LineWidth',linesize);
    %
    %     plot3(AgentNom(jj).X_LyapTrack(1,1),...
    %         AgentNom(jj).X_LyapTrack(1,2),AgentNom(jj).X_LyapTrack(1,3),...
    %         's','LineWidth',4,'MarkerEdgeColor',cMat0(ij,:),...
    %         'MarkerFaceColor',cMat0(ij,:),'MarkerSize',15);
    %
    %     plot3(AgentNom(jj).X_LyapTrack(end,1),...
    %         AgentNom(jj).X_LyapTrack(end,2),AgentNom(jj).X_LyapTrack(end,3),...
    %         'p','LineWidth',4,'MarkerEdgeColor',cMat0(ij,:),...
    %         'MarkerFaceColor',cMat0(ij,:),'MarkerSize',20);
    %
    %      arrow3d(AgentNom(jj).ArrowsBase,AgentNom(jj).ArrowsHead,[c6;c6],[0.65,0.85]);
    %
    % end
    %
    %
    % aa = textcase(1);
    % for jj = 1:length(AgentNom(aa).IterationIndex)
    %     indxc = AgentNom(aa).IterationIndex(jj);
    %     vel = AgentNom(textcase(1)).DesiredState(indxc,4:6)'; b3 = vel/norm(vel); zhat = [0;0;1];
    %     b2 = cross(b3,zhat); b2 = b2/norm(b2);
    %     b1 = cross(b2,b3);   b1 = b1/norm(b1);
    %     R = [b1 b2 b3];
    %     plt2 = MotionCone(AgentNom(aa).KeepOutRadius(jj),ConeAngle,R,AgentNom(aa).DesiredState(indxc,1:3),c11);
    %
    %     if jj == 1
    %         plot3(AgentNom(5).DesiredState(indxc,1),...
    %         AgentNom(5).DesiredState(indxc,2),...
    %         AgentNom(5).DesiredState(indxc,3),'bo',...
    %         'LineWidth',4,...
    %         'MarkerEdgeColor',c12,...
    %         'MarkerFaceColor',c12,...
    %         'MarkerSize',10);
    %     else
    %     plt5 = plot3(AgentNom(6).DesiredState(indxc,1),...
    %         AgentNom(6).DesiredState(indxc,2),...
    %         AgentNom(6).DesiredState(indxc,3),'bo',...
    %         'LineWidth',4,...
    %         'MarkerEdgeColor',c12,...
    %         'MarkerFaceColor',c12,...
    %         'MarkerSize',10);
    %     end
    % end
    %
    % hL = legend([plt1(end,3), plt1(end,2), plt1(end,1), plt0,plt5,plt2],...
    %     {"Agent 6's Path","Agent 5's Path"...
    %     'Agent 4 (Initial)',...
    %     'Agent 4 (Updated)',...
    %     'Detected Obstacles',...
    %     'Avoidance Cone'},...
    %     'AutoUpdate','off');
    % % saveLegendToImage(fh2, hL, './Figures/3DLegendAgent4')
    % hL.FontSize = 30;
    % newPosition = [0.72 0.35 0.1 0.1];
    % newUnits = 'normalized';
    % set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1);
    %
    % set(gca,'FontSize',40)
    % % set all units inside figure to normalized so that everything is scaling accordingly
    % set(findall(fh2,'Units','pixels'),'Units','normalized');
    % % do show figure on screen
    % set(fh2, 'visible', 'on')
    % % set figure units to pixels & adjust figure size
    % fh2.Units = 'pixels';
    % fh2.OuterPosition = [10 10 1000 800]; % [10 10 1400 1000];
    %
    % % res = 500; % resolution of figure to be saved in dpi
    % % set(fh2,'PaperPositionMode','manual')
    % % fh2.PaperUnits = 'inches';
    % % fh2.PaperPosition = [0 0 8580 6320]/res;
    % % print(fh2,'./Figures/3DPlotAgent4','-dpng',sprintf('-r%d',res))
end
% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% ######################## Required Functions #############################
% -------------------------------------------------------------------------
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   One Agent   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% --------------------------- PDE for pdepe ------------------------------
function [c,f,s] = mypdexpde(x,t,u,DuDx,k,N) % Define PDE; right-hand-side of AGHF
XChief = full(ChiefState(x));
LL = full(L(u,DuDx,k,x,XChief));
CasADiresult = full(dLdx(u,DuDx,k,x,XChief,LL))';
NN = length(u);
CasADir = [CasADiresult(1:NN), CasADiresult(NN+1:2*NN)];

f =  CasADir(:,2); 
s = -CasADir(:,1); 
c = ones(N,1);

end

% --------------------- Initial condition for pdepe ----------------------
function u0 = mypdexic(x, X0, Xf, T)    % Initial condition for AGHF
% Sinusoidal initial condition
freq1 = pi/2 + 2*pi;
freq2 = pi + 2*pi;

u0 = X0.*cos(freq1.*x./T) ...
    + Xf.*sin(freq1.*x./T) ...
    - (X0+Xf).*sin(freq2.*x./T);

end

% --------------------- Boundary condition for pdepe ---------------------
function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t, X0, Xf, N) % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!   Multiple Agents   !!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% --------------------------- PDE for pdepe ------------------------------
function [c,f,s] = ReplanningPDE(x,t,u,DuDx,k,N,ObstaclesList,kb,pb,AgentNom)

Nobs      = length(ObstaclesList);
XChief    = full(ChiefState(x));
Obstacles = zeros(3,9);
for j = 1:Nobs
    ObsID = ObstaclesList(j);
    ObstacleStates = AgentNom(ObsID).ObstaclePredictionInterpolant(x)';
    Obstacles(:,j) = ObstacleStates(1:3,:);
end

[Lr1,Lr2,Lr3,Lr4,Lr5,Lr6,Lr7,Lr8,Lr9] = Lreplan(u,DuDx,k,x,XChief,Obstacles,kb,pb);
DL = Lreplandx(u,DuDx,k,x,XChief,Obstacles,kb,pb,Lr1,Lr2,Lr3,Lr4,Lr5,Lr6,Lr7,Lr8,Lr9);

if Nobs == 1
    CasADir = [full(DL(1,1:6))', full(DL(1,7:12))'];
elseif Nobs == 2
    CasADir = [full(DL(2,1:6))', full(DL(2,7:12))'];
elseif Nobs == 3
    CasADir = [full(DL(3,1:6))', full(DL(3,7:12))'];
elseif Nobs == 4
    CasADir = [full(DL(4,1:6))', full(DL(4,7:12))'];
elseif Nobs == 5
    CasADir = [full(DL(5,1:6))', full(DL(5,7:12))'];
elseif Nobs == 6
    CasADir = [full(DL(6,1:6))', full(DL(6,7:12))'];
elseif Nobs == 7
    CasADir = [full(DL(7,1:6))', full(DL(7,7:12))'];
elseif Nobs == 8
    CasADir = [full(DL(8,1:6))', full(DL(8,7:12))'];
elseif Nobs == 9
    CasADir = [full(DL(9,1:6))', full(DL(9,7:12))'];
end
    
% dL/dx
pLx  = CasADir(:,1); % EL(:,1); %
% dL/d(dot_x)
pLxd = CasADir(:,2); % EL(:,2); %

f = pLxd;
s = -pLx;
c = ones(N,1);

end % Define PDE; right-hand-side of AGHF


% --------------------- Initial condition for pdepe ----------------------
function u0 = ReplanningIC(x,AgentNom,AgentID,To,X0,Xf,T) % Initial condition for AGHF  X0, Xf, To, T,

u0 = AgentNom(AgentID).DesiredStateInterpolant(x)';

end

% --------------------- Boundary condition for pdepe ---------------------
function [pl,ql,pr,qr] = ReplanningBCs(xl,ul,xr,ur,t, X0, Xf, N) % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!   Nondimetionalized Dynamics Models   !!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ------------------- Nondimetionalized Dynamics Models -------------------
function [Xdot] = NonDimentionalized_NLode(t,X)
global mu_Earth Tc
Xchief = full(ChiefState(t));
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
rt = Xchief(1:3);  vt = Xchief(4:6); rt_norm = (rt.'*rt).^(1/2);
% Redimensionalizing the variables in the ODE
n = 6; m = length(X)/n;
Xnew = reshape(X,[n,m]);
RhoMat = Xnew(1:3,:); RhoDotMat = Xnew(4:6,:);
RcMat = [rt_norm; 0; 0] + RhoMat;
RcMatNorm = vecnorm(RcMat,2,1);

% Calculating the respective accelerations
TN = DCM(rt,vt); % DCM's between target's RIC frame and ECI frame

% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(rt)*vt;
Omega = TN*( h_vec/rt_norm^2);
Omegadot = TN*( -2*(rt.'*vt)*h_vec/rt_norm^4 );

% relative gravitationall acceleration and coriolis effects
DelAgMat     = -mu_Earth./(RcMatNorm.^3).*RcMat + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
RotEffectMat = - 2*tilde(Omega)*RhoDotMat ...
               - tilde(Omega)*(tilde(Omega)*RhoMat) ...
               - tilde(Omegadot)*RhoMat;
% Nondimetional relative ODE of the deputy
XdotMat = [RhoDotMat; (DelAgMat + RotEffectMat)];
Xdot = XdotMat(:)*Tc;

end
function [Xdot] = NonDimentionalized_NLode0(t,X)
global mu_Earth Tc
Xchief = full(ChiefState0(t));
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
rt = Xchief(1:3);  vt = Xchief(4:6); rt_norm = (rt.'*rt).^(1/2);
% Redimensionalizing the variables in the ODE
n = 6; m = length(X)/n;
Xnew = reshape(X,[n,m]);
RhoMat = Xnew(1:3,:); RhoDotMat = Xnew(4:6,:);
RcMat = [rt_norm; 0; 0] + RhoMat;
RcMatNorm = vecnorm(RcMat,2,1);

% Calculating the respective accelerations
TN = DCM(rt,vt); % DCM's between target's RIC frame and ECI frame

% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(rt)*vt;
Omega = TN*( h_vec/rt_norm^2);
Omegadot = TN*( -2*(rt.'*vt)*h_vec/rt_norm^4 );

% relative gravitationall acceleration and coriolis effects
DelAgMat     = -mu_Earth./(RcMatNorm.^3).*RcMat + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
RotEffectMat = - 2*tilde(Omega)*RhoDotMat ...
               - tilde(Omega)*(tilde(Omega)*RhoMat) ...
               - tilde(Omegadot)*RhoMat;
% Nondimetional relative ODE of the deputy
XdotMat = [RhoDotMat; (DelAgMat + RotEffectMat)];
Xdot = XdotMat(:)*Tc;
end
function [Xdot] = NonDimentionalized_NLodef(t,X)
global mu_Earth Tc
Xchief = full(ChiefStatef(t));
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
rt = Xchief(1:3);  vt = Xchief(4:6); rt_norm = (rt.'*rt).^(1/2);
% Redimensionalizing the variables in the ODE
n = 6; m = length(X)/n;
Xnew = reshape(X,[n,m]);
RhoMat = Xnew(1:3,:); RhoDotMat = Xnew(4:6,:);
RcMat = [rt_norm; 0; 0] + RhoMat;
RcMatNorm = vecnorm(RcMat,2,1);

% Calculating the respective accelerations
TN = DCM(rt,vt); % DCM's between target's RIC frame and ECI frame

% Computing the angular velocity and angular acceleration of the target frame
h_vec = tilde(rt)*vt;
Omega = TN*( h_vec/rt_norm^2);
Omegadot = TN*( -2*(rt.'*vt)*h_vec/rt_norm^4 );

% relative gravitationall acceleration and coriolis effects
DelAgMat     = -mu_Earth./(RcMatNorm.^3).*RcMat + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
RotEffectMat = - 2*tilde(Omega)*RhoDotMat ...
               - tilde(Omega)*(tilde(Omega)*RhoMat) ...
               - tilde(Omegadot)*RhoMat;
% Nondimetional relative ODE of the deputy
XdotMat = [RhoDotMat; (DelAgMat + RotEffectMat)];
Xdot = XdotMat(:)*Tc;

end
function [Xdot] = NonDimentionalized_NLodeCasADi(t,X,XChief)
global mu_Earth Tc
tilde = @(v) [0    -v(3)  v(2);
              v(3)   0   -v(1);
              -v(2)  v(1)    0];
% Collecting relevant quantities
rt = XChief(1:3); vt = XChief(4:6); rt_norm = (rt.'*rt).^(1/2);
% Redimensionalizing the variables in the ODE
n = 6; m = length(X)/n; DelAgMat = [];
Xnew   = reshape(X,[n,m]);
RhoMat = Xnew(1:3,:); RhoDotMat = Xnew(4:6,:);
RcMat  = [rt_norm; 0; 0] + RhoMat;

% Calculating the respective accelerations
TN = DCM(rt,vt); % DCM's between target's RIC frame and ECI frame

% Computing the angular velocity and angular acceleration of the target frame
h_vec    = tilde(rt)*vt;
Omega    = TN*(h_vec/rt_norm^2);
Omegadot = TN*(-2*(rt.'*vt)*h_vec/rt_norm^4);

% relative gravitationall acceleration and coriolis effects
for i = 1:m
    DelAgMat = [DelAgMat, -mu_Earth/(RcMat(:,i).'*RcMat(:,i)).^(3/2)*RcMat(:,i)...
                + mu_Earth/rt_norm^2*[1 0 0]']; % (2-body) gravitatinal];
end
RotEffectMat = - 2*tilde(Omega)*RhoDotMat ...
               - tilde(Omega)*(tilde(Omega)*RhoMat) ...
               - tilde(Omegadot)*RhoMat;
% Nondimetional relative ODE of the deputy
XdotMat = [RhoDotMat; (DelAgMat + RotEffectMat)];
Xdot    = XdotMat(:)*Tc;

end









