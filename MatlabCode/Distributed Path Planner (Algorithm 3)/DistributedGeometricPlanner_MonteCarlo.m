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
global mu_Earth tspan mu_dot ...
    Num_Agent a_chief e_chief Tc
% ----------------- Calculating the boundaries conditions ----------------
mu_Earth   = 3.986004415E5; % Earth gravitational parameter
Num_Agent  = 6;             % Number of agents in the system
N          = 6;             % Size of the state space
McRuns     = 1;           % Number of Monte Calo simulations 
MonteCarlo = struct('CrtStd',cell(1),...
                    'NumbReplanningAgent',cell(1),...
                    'ReplanningSequences',cell(1));

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
Ttraj0     = linspace(0,Period,1e4);
Ttrajf     = linspace(IntTime,Period+IntTime,1e4);
options = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-17);
mu_dot  = sqrt(mu_Earth/a_chief^3); % Chief's mean motion
Rc  = a_chief*eye(3); Tc = Period;
inc_chief = deg2rad(inc_chief); BigOmg_chief = deg2rad(BigOmg_chief);
LitOmg_chief = deg2rad(LitOmg_chief); M_chief1 = deg2rad(M_chief1);
COE1 = [a_chief,e_chief,inc_chief,BigOmg_chief,LitOmg_chief,M_chief1];
[Position_target1,Velocity_target1] = COEstoRV(COE1,mu_Earth);
X0_Chief1 = [Position_target1; Velocity_target1];
[~, Xnom] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),tspan,X0_Chief1,options);
[~, Xnom0] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),Ttraj0,Xnom(1,:)',options);
[~, Xnomf] = ode113(@(t,X)M2BodyOde(t,X,mu_Earth),Ttrajf,Xnom(end,:)',options);
ChiefMeanMotion = sqrt(mu_Earth/a_chief^3);


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!   Monte Carlo Simulation   !!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for Mc = 1:McRuns
    % ======================================================================= %
    % --------------- Setting the deputy initial conditions ------------------
    % ======================================================================= %
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
                
                dele1     = -1/(3*a_chief) *randi([1 9]);
                deli1     = -((-1)^k)/(a_chief) *randi([1 9]);
                delLitOmg = (1/5)*pi*1E-5 *randi([1 9]);
                delBigOmg = (1/5)*pi*1E-5 *randi([1 9]);
                delM      = (1/5)*pi*1E-6 *randi([1 9]);
                delq1_1   = dele1*cos(LitOmg_chief) - e_chief*sin(LitOmg_chief)*delLitOmg;
                delq2_1   = dele1*sin(LitOmg_chief) + e_chief*cos(LitOmg_chief)*delLitOmg;
                deltheta  = delLitOmg + delM; % Relative true latitude in rad
                delCOE_0  = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
                X0(idx,1) = AMap1*delCOE_0;
                
                dele1     = 1/(2*a_chief) *randi([1 9]);
                deli1     = ((-1)^k)/(a_chief) *randi([1 9]);
                delLitOmg = (1/5)*pi*1E-5 *randi([1 9]);
                delBigOmg = (1/5)*pi*1E-5 *randi([1 9]);
                delM      = (1/5)*pi*1E-6 *randi([1 9]);
                delq1_1   = dele1*cos(LitOmg_chief) - e_chief*sin(LitOmg_chief)*delLitOmg;
                delq2_1   = dele1*sin(LitOmg_chief) + e_chief*cos(LitOmg_chief)*delLitOmg;
                deltheta  = delLitOmg + delM; % Relative true latitude in rad
                delCOE_f  = [dela, deltheta, deli1, delq1_1, delq2_1, delBigOmg]';
                Xf(idx,1) = AMap2*delCOE_f;
            end
            
            X02 = X0; shuffle1 = randperm(N,N);
            for k = 1:length(shuffle1)
                j = shuffle1(k); idx  = 1+N*(k-1):N*k; idx2 = 1+N*(j-1):N*j;
                X0(idx,1) = X02(idx2,1);
            end
            Xf2 = Xf; shuffle2 = randperm(N,N);
            for k = 1:length(shuffle2)
                j = shuffle2(k); idx  = 1+N*(k-1):N*k; idx2 = 1+N*(j-1):N*j;
                Xf(idx,1) = Xf2(idx2,1);
            end
        end % Initial conditions to the Monte Carlo simulation
        
        [~, Traj1] = ode113(@(t,x)NonDimentionalized_NLode0(t,x),TT0,X0,options);
        [~, Traj2] = ode113(@(t,x)NonDimentionalized_NLodef(t,x),TTf,Xf,options);
        AgentNom = struct('Xrel0',{Traj1(:,1:6),  Traj1(:,7:12),...
            Traj1(:,13:18),Traj1(:,19:24),...
            Traj1(:,25:30),Traj1(:,31:36)},...
            'Xrelf',{Traj2(:,1:6),  Traj2(:,7:12),...
            Traj2(:,13:18),Traj2(:,19:24),...
            Traj2(:,25:30),Traj2(:,31:36)});
    end
    
    
    % ======================================================================= %
    % --------------------------- AGHF parameters ----------------------------
    % ======================================================================= %
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
    
    
    % ======================================================================= %
    % ------- Solving the Unconstraint Problem (Assuming No Obstacles) -------
    % ======================================================================= %
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
    
    
    % ======================================================================= %
    % -- Tracking the Computed Optimal Path (Lyapunov Tracking Controller) ---
    % ======================================================================= %
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
    
    % ======================================================================= %
    % -------------------- Receding Horizon Path Planner ---------------------
    % ======================================================================= %
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
                        AgentNom(ObstaclesList(ijk)).ObstaclePredictionInterpolant = ...
                            griddedInterpolant(xnew(idx:end),AgentNom(ijk).DesiredState(idx:end,:),'spline');
                    end
                    
                    %=========================== Replanning sequence ===========================%
                    AgentX0 = AgentNom(j).DesiredState(idx,:)'; AgentXf = Xf(1+N*(j-1):N*j);
                    pb(1:length(ObstaclesList),:) = [Separation(ObstaclesList,j)/2 ObstaclesList'];
                    kb(1:length(ObstaclesList),:) = [kbMat(ObstaclesList)' ObstaclesList'];
                    tic
                    NewPath = pdepe(m,@(x,t,u,DuDx) ReplanningPDE(x,t,u,DuDx,Penalty,N,ObstaclesList,kb,pb,AgentNom),...
                        @(x) ReplanningIC(x,AgentNom,j),...
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
        
        ijk = 0;
        for k = 1:Num_Agent
            AgentNom(k).TotalCtrl_Int = 1e6*simpsons(vecnorm(AgentNom(k).Crtl,2,2),0,0.75,[]);
            AgentNom(k).TotalCtrl_New = 1e6*simpsons(vecnorm(AgentNom(k).CrtlNew,2,2),0,0.75,[]);
            MonteCarlo.CrtStd(Mc,k) = 100*(AgentNom(k).TotalCtrl_New - AgentNom(k).TotalCtrl_Int)/AgentNom(k).TotalCtrl_Int;
            MonteCarlo.ReplanningSequences(Mc,k) = length(AgentNom(k).IterationIndex);
            if ~isempty(AgentNom(k).AGHFCostNew)
                ijk = ijk+1;
            end
        end
        MonteCarlo.NumbReplanningAgent = ijk;
    end
    
end
save MonteCarloResults.mat MonteCarlo -v7.3 % Saving the Results

number = '2028171774';
carrier = 'AT&T';
message = 'Your MatLab code is done running!';
send_text_message(number,carrier,message)











