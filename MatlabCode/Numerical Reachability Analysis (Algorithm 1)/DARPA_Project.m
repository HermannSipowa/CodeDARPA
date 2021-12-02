%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Hermann Kaptui Sipowa %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
clear
close all
clc
start_up
format long e
clc
global JD mu_Earth AU Target ae
ae       = 6378.136;
mu_Earth = 3.986004415E5;
JD       = 2456296.25;
AU       = 149597870.7;
Tsize    = 1e2;
SimTime  = 1;
Target   = Spacecraft3D([5, 5, 1, .25, .5, .5, 0.2, 0.5, 1e-2, 2, 1.5]); % Target characteristic
options  = odeset('RelTol',2.22045e-14,'AbsTol',2.22045e-30);

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
c12 = rgb('Black');

addpath('../casadi-osx-matlabR2015a-v3.5.5')
import casadi.*
CasADiopts = struct('main', true, 'mex', true);

%% =========== Setting the inital condition of the simulation ========== %%
anom   = ae+650; %4.2164e4;    % Semi-major axis in Km  % ae+600;
e      = 0.02;        % Eccentricity
inc    = deg2rad(50); % Inclination in deg0
BigOmg = deg2rad(10); % RAAN in deg
LitOmg = deg2rad(10); % AOP in deg
M      = deg2rad(0);  % Initial mean anomaly

Period  = (2*pi)*sqrt(anom^3/mu_Earth);
IntTime = 2*Period;
tspan   = linspace(0,IntTime,Tsize);

COE    = [anom,e,inc,BigOmg,LitOmg,M];
[Position_target1,Velocity_target1] = COEstoRV(COE,mu_Earth);
qr_chasser0 = [8 3 7 9]'; qr_chasser0 =  qr_chasser0/norm(qr_chasser0);
Omega_chaser0 = 2e-2*[deg2rad(-.5) deg2rad(.4) -deg2rad(.35)]'; % Angular velocity in inertial frame
CN = Quaternion_to_DCM(qr_chasser0); Omega_chaser_B0 = CN*Omega_chaser0;
Xnom0 = [Position_target1; Velocity_target1];
Xsrp1 = [qr_chasser0; Omega_chaser0;...
    Position_target1; Velocity_target1;...
    qr_chasser0; Omega_chaser0];
Xsrp2 = [qr_chasser0; Omega_chaser0;...
    Position_target1; Velocity_target1];

%% ------------------------ AutoDiff using CasADi -------------------------
for l = 1
    r = SX.sym('r',3); v = SX.sym('v',3);
    q = SX.sym('q',4); t = SX.sym('t',1);
    % Switching condition for velocity tracker
    qt = q/norm(q);
    Udrag1 = Drag(r,v,qt,Target);
    [u_srp, ~, ~] = FlatPlate_srp(t, r, Target, qt);
    BN = Quaternion_to_DCM(qt);
    Usrp1 = BN\u_srp;
    LVLH = DCM(r,v); Fdist = LVLH*(Usrp1+u_srp);
    OBE = RVtoCOEs_MeanAnomaly(r,v,mu_Earth);
    a = OBE(1); e = OBE(2); f = OBE(6);
    h = norm(cross(r,v)); p = a*(1-e^2);
    
    da_dt = Function('da_dt',{r,v,q,t},...
        {2*a^2/h * [e*sin(f), p/norm(r), 0] * Fdist},...
        {'r','v','q','t',},...
        {'adot'});
    dda_dt = da_dt.jacobian();
    
    da_dt.generate('da_dt.c',CasADiopts);
    dda_dt.generate('dda_dt.c',CasADiopts);
    
    % Generating Mex files from the c files
    mex da_dt.c -largeArrayDims
    mex dda_dt.c -largeArrayDims
    clear da_dt dda_dt fdrift
end

%%
[~, Xsrp] = ode113(@(t,X)SRPandDragODE(t,X,Target),tspan,Xsrp1);
[~, Xnom] = ode113(@(t,X)M2BPODE(t,X),tspan,Xnom0,options);
[~, Xsrpony] = ode113(@(t,X)SRPonlyODE(t,X,Target),tspan,Xsrp2);

XDiff = Xnom - Xsrp(:,8:13);
XDiffSRPonly = Xnom - Xsrpony(:,8:13);



%%
for ll = 1
    a_srp = nan(Tsize,1);
    a_srponly = nan(Tsize,1); dadtonly = nan(Tsize,1); Nu = nan(Tsize,1);
    dadt = nan(Tsize,1); % addot = nan(Tsize,1);
    for i = 1:Tsize
        t = tspan(i);
        % ======================================== %
        r = Xsrp(i,8:10)';
        v = Xsrp(i,11:13)';
        OBE = RVtoCOEs_MeanAnomaly(r,v,mu_Earth);
        a_srp(i) = OBE(1) - anom;
        e_srp(i) = OBE(2);
        % =================================================================== %
        q = Xsrp(i,1:4)';
        qt = q/norm(q);
        Udrag1 = Drag(r,v,qt,Target);
        [u_srp, ~, ~] = FlatPlate_srp(t, r, Target, qt);
        BN = Quaternion_to_DCM(qt);
        Usrp1 = BN\u_srp;
        LVLH = DCM(r,v); Fdist = LVLH*(Usrp1+Udrag1);
        OBE = RVtoCOEs_MeanAnomaly(r,v,mu_Earth);
        a = OBE(1); e = OBE(2); f = OBE(6);
        h = norm(cross(r,v)); p = a*(1-e^2);
        dadt(i) = 2*a^2/h * [e*sin(f), p/norm(r), 0] * Fdist;
        % ad = dda_dt(r,v,q,t,dadt(i)); addot(i) = ad(end);
    
        % =================================================================== %
        r = Xsrpony(i,8:10)';
        v = Xsrpony(i,11:13)';
        OBE = RVtoCOEs_MeanAnomaly(r,v,mu_Earth);
        a_srponly(i) = OBE(1) - anom;
    
        % Computing SRP acceleration (Flat Plate model)
        q = Xsrpony(i,1:4)'/norm(Xsrpony(i,1:4)');
        [u_srp, ~, ~] = FlatPlate_srp(t, r, Target, q);
        BN = Quaternion_to_DCM(q);
        Usrp1 = BN\u_srp;
        Udrag1 = Drag(r,v,q,Target);
        % Switching condition for the controller
        LVLH = DCM(r,v); Fdist = LVLH*Usrp1;
        OBE = RVtoCOEs_MeanAnomaly(r,v,mu_Earth);
        a = OBE(1); e = OBE(2); f = OBE(6);
        h = norm(cross(r,v)); p = a*(1-e^2);
        dadtonly(i) = 2*a^2/h * [e*sin(f), p/norm(r), 0] * Fdist;
        Nu(i) = f;
    
    end
    
   
    % ====================== Plotting the results ====================== %
    close all
    Time = tspan/3600;
    Nu   = unwrap(Nu);
    figure
    plot(Time,e_srp,'Color',c6,'LineWidth',3);
    grid on
    grid minor
    
    
    fh2 = figure;
    subplot(2,1,1)
    plt1 = plot(Time,a_srp,'Color',c6,'LineWidth',3);
    hold on
    plt2 = plot(Time,a_srponly,'Color',c5,'LineWidth',3);
    yl = yline(0,'-.','h \approx 500 km','LineWidth',3);
    ylim([min(a_srp)-(max(a_srp))/5, inf])
    xlim([min(Time), max(Time)+0.01])
    yl.LabelHorizontalAlignment = 'right';
    yl.LabelVerticalAlignment = 'top';
    yl.Color = c7;
    yl.FontSize = 20;
    yl.FontWeight = 'bold';
    xlabel('time [hr]')
    grid on
    grid minor
    ylabel('$\delta a$ [km]')
    set(gca,'FontSize',30)
    title('semi-major axis change')
    
    subplot(2,1,2)
    plot(Nu,dadt,'Color',c6,'LineWidth',3)
    hold on
    plot(Nu,dadtonly,'Color',c5,'LineWidth',3)
    yl = yline(0,'-.','LineWidth',3);
    ylim([min(min(dadtonly),min(dadt))-(min(max(dadtonly),max(dadt)))/5,...
        max(max(dadtonly),max(dadt))+(min(max(dadtonly),max(dadt)))/5])
    xlim([min(Nu), max(Nu)+.01])
    yl.Color = c7;
    yl.FontSize = 10;
    yl.FontWeight = 'bold';
    xlabel('$f$ [rad]')
    ylabel('$\frac{da}{dt}$[km/s]','interpreter','latex')
    pmin = min(Nu(:));
    pmax = max(Nu(:));
    pimin = floor(pmin/pi);
    pimax = ceil(pmax/pi);
    xticks((pimin:1/2:pimax) * pi);
    xticklabels( {'0','$\frac{1}{2}\pi$','$\pi$','$\frac{3}{2}\pi$',...
        '$2\pi$','$\frac{5}{2}\pi$','$3\pi$','$\frac{7}{2}\pi$','$4\pi$'} )
    set(gca,'FontSize',30)
    grid on
    grid minor
    
    hL = legend([plt1, plt2],...
        {'SRP + Drag','SRP only'},'AutoUpdate','off');
    hL.FontSize = 20;
    newPosition = [0.25 0.8 0.05 0.05];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1  );
    
    set(findall(fh2,'Units','pixels'),'Units','normalized');
    % do show figure on screen
    set(fh2, 'visible', 'on')
    % set figure units to pixels & adjust figure size
    fh2.Units = 'pixels';
    fh2.OuterPosition = [10 10 1000 800];
    
    res = 500;
    % recalculate figure size to be saved
    set(fh2,'PaperPositionMode','manual')
    fh2.PaperUnits = 'inches';
    fh2.PaperPosition = [0 0 6580 5320]/res;
    % save figure
    % print(fh2,'SRPandDrag','-dpng',sprintf('-r%d',res))
end


%% ============ Monte Cacrlo =========== %%
% =========== Setting the inital condition of the simulation ========== %%
McRuns  = 50;
a = ae + 233 + 50*(0:1:McRuns);
e      = 0.02;        % Eccentricity
inc    = deg2rad(50); % Inclination in deg0
BigOmg = deg2rad(10); % RAAN in deg
LitOmg = deg2rad(10); % AOP in deg
M      = deg2rad(0);  % Initial mean anomaly
Target   = Spacecraft3D([5, 5, 1, .25, 2.5, .5, 0.2, 0.5, 1e-2, 2 1.5]); % Target characteristic
Period  = (2*pi)*sqrt(a(end)^3/mu_Earth);
IntTime = 2*Period;
tspan   = linspace(0,IntTime,Tsize);
qr_chasser0 = [8 3 7 9]'; qr_chasser0 =  qr_chasser0/norm(qr_chasser0);
Omega_chaser0 = 2e-2*[deg2rad(-.5) deg2rad(.4) -deg2rad(.35)]'; % Angular velocity in inertial frame
a_srp = nan(Tsize,McRuns);

% Monte Carlo on the different altitudes
for j = 1:McRuns
    COE  = [a(j),e,inc,BigOmg,LitOmg,M];
    [Position_target1,Velocity_target1] = COEstoRV(COE,mu_Earth);
    Xsrp1 = [qr_chasser0; Omega_chaser0;...
             Position_target1; Velocity_target1;...
             qr_chasser0; Omega_chaser0];
    [~, Xsrp] = ode113(@(t,X)SRPandDragODE(t,X,Target),tspan,Xsrp1);
    h(j) = a(j)*(1-e) - ae; % Altitude at periasis
    for i = 1:Tsize
        r = Xsrp(i,8:10)';
        v = Xsrp(i,11:13)';
        OBE = RVtoCOEs_MeanAnomaly(r,v,mu_Earth);
        a_srp(i,j) = OBE(1) - a(j);
    end
end
% save Asrp10.mat a_srp -v7.3

%%
a_srp = load('Asrp10.mat'); a_srp = a_srp.a_srp;
% Plotting the results
close all
Time = tspan/3600;
plt  = zeros(McRuns);
jetcustom = jet(McRuns);
fh2 = figure;
for j = 7:McRuns
    plt(:,j)  = plot(Time,a_srp(:,j),'LineWidth',3, 'Color',  jetcustom(j,:));
    hold on
end
yl = yline(0,'-.','Initial semi-major axis','LineWidth',3);
xlim([min(Time), max(Time)+0.01])
yl.LabelHorizontalAlignment = 'right';
yl.LabelVerticalAlignment = 'top';
yl.Color = c7;
yl.FontSize = 30;
yl.FontWeight = 'bold';
xlabel('time [hr]')
grid on
grid minor
ylabel('$\delta a$ [km]')

colormap(jet(McRuns))
cb = colorbar();
caxis([h(1) h(end)])
ylabel(cb,'Altitude [km]','FontSize',40)

set(gca,'FontSize',30)
title('semi-major axis change')

% hL = legend([plt1, plt2],...
%     {'SRP + Drag','SRP only'},'AutoUpdate','off');
% hL.FontSize = 20;
% newPosition = [0.25 0.75 0.05 0.05];
% newUnits = 'normalized';
% set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1  );

set(findall(fh2,'Units','pixels'),'Units','normalized');
% do show figure on screen
set(fh2, 'visible', 'on')
% set figure units to pixels & adjust figure size
fh2.Units = 'pixels';
fh2.OuterPosition = [10 10 1200 900];

res = 500;
% recalculate figure size to be saved
set(fh2,'PaperPositionMode','manual')
fh2.PaperUnits = 'inches';
fh2.PaperPosition = [0 0 8580 5320]/res;
% save figure
% print(fh2,'SRPandDragMonteCarlo25','-dpng',sprintf('-r%d',res))

%% Monte carlo on different area to mass ratios
McRuns = 25;
a = ae + 633;
mass = 1:1:McRuns;
a_AMR = nan(Tsize,McRuns);
for j = 1:McRuns
    j
    Target   = Spacecraft3D([5, 5, 1, .25,...
               mass(j), 0, 0.2, 0.5, 1e-2, 2 1.5]); % Target characteristic
    COE  = [a,e,inc,BigOmg,LitOmg,M];
    [Position_target1,Velocity_target1] = COEstoRV(COE,mu_Earth);
    Xsrp1 = [qr_chasser0; Omega_chaser0;...
             Position_target1; Velocity_target1;...
             qr_chasser0; Omega_chaser0];
    [~, Xsrp] = ode113(@(t,X)SRPandDragODE(t,X,Target),tspan,Xsrp1);
    
    for i = 1:Tsize 
        r = Xsrp(i,8:10)';
        v = Xsrp(i,11:13)';
        OBE = RVtoCOEs_MeanAnomaly(r,v,mu_Earth);
        a_AMR(i,j) = OBE(1) - a;
    end
end
 save AreaToMassRatio.mat a_AMR -v7.3

close all
cMat = [c1;c5;c6;c11];
cMat2 = [c2;c3;c5;c9];
Time = tspan/3600;
plt  = zeros(McRuns);
jetcustom = jet(McRuns);
fh2 = figure;
for j = 1:McRuns
    plt(:,j)  = plot(Time,a_AMR(:,j),'LineWidth',3, 'Color',  jetcustom(j,:));
    hold on
end
yl = yline(0,'-.','h \approx 400 km','LineWidth',3);
xlim([min(Time), max(Time)+0.01])
yl.LabelHorizontalAlignment = 'right';
yl.LabelVerticalAlignment = 'top';
yl.Color = c7;
yl.FontSize = 30;
yl.FontWeight = 'bold';
xlabel('time [hr]')
grid on
grid minor
ylabel('$\delta a$[km]')

colormap(jet(McRuns))
cb = colorbar();
caxis([mass(end) mass(end)])
oldcmap = colormap;
colormap( flipud(oldcmap) );
ylabel(cb,'Araes/Mass','FontSize',40)

set(gca,'FontSize',30)
title('semi-major axis change')

% hL = legend([plt1, plt2],...
%     {'SRP + Drag','SRP only'},'AutoUpdate','off');
% hL.FontSize = 20;
% newPosition = [0.25 0.75 0.05 0.05];
% newUnits = 'normalized';
% set(hL,'Position', newPosition,'Units', newUnits,'NumColumns',1  );

set(findall(fh2,'Units','pixels'),'Units','normalized');
% do show figure on screen
set(fh2, 'visible', 'on')
% set figure units to pixels & adjust figure size
fh2.Units = 'pixels';
fh2.OuterPosition = [10 10 1200 900];

res = 500;
% recalculate figure size to be saved
set(fh2,'PaperPositionMode','manual')
fh2.PaperUnits = 'inches';
fh2.PaperPosition = [0 0 8580 5320]/res;
% save figure
% print(fh2,'AeraMassRatioMonteCarlo','-dpng',sprintf('-r%d',res))

%% -------------------------------------------------------------------------
% ######################## Required Functions #############################
% -------------------------------------------------------------------------
function [Xdot] = M2BPODE(t,X)
global mu_Earth
Xdot(1:3,1)  = X(4:6,1);
Xdot(4:6,1) = -mu_Earth/norm(X(1:3,:))^3 * X(1:3,1);
end
% -------------------------------------------------------------------------
function Xdot = SRPonlyODE(t, X, Target)
global mu_Earth
q = X(1:4)/norm(X(1:4));
Omega = X(5:7);
r = X(8:10); v = X(11:13);
r_norm = norm(r);
tilde = @(v) [0     -v(3)  v(2) ;
              v(3)   0    -v(1) ;
             -v(2)   v(1)  0   ];
Moment_Inertia = diag(Target.Moment_Of_Inertia_Calculator());

% Rotational state differential equation
W_Matrix = [0, -Omega(1), -Omega(2), -Omega(3);
    Omega(1),  0,  Omega(3), -Omega(2);
    Omega(2), -Omega(3),  0,  Omega(1);
    Omega(3), Omega(2), -Omega(1), 0];


% Computing the forcing accelerations acting on the spacecraft
% Udrag = Drag(r,v,q,Target);
[u_srp, XS, VS] = FlatPlate_srp(t, r, Target, q);
BN = Quaternion_to_DCM(q);
Usrp1 = BN\u_srp;

% Switching condition for the controller
LVLH = DCM(r,v); Fdist = LVLH*(Usrp1);
OBE = RVtoCOEs_MeanAnomaly(r,v,mu_Earth);
a = OBE(1); e = OBE(2); f = OBE(6);
h = norm(cross(r,v)); p = a*(1-e^2);
da_dt = 2*a^2/h * [e*sin(f), p/norm(r), 0] * Fdist;

% Compting the torque required to track the sun
P = 30*diag([1 2 3]); k = 30; % Controller gains
U_torque = u_Cframe_QuaternionTracking(P,k,Moment_Inertia,XS,VS,r,v,q,Omega);

if da_dt<0
    Usrp = zeros(3,1);
else
    Usrp = Usrp1;
end

% Integrating the the target trajectory
Xdot(1:4) = 1/2*W_Matrix*q;
Xdot(5:7) = Moment_Inertia\(-tilde(Omega)*(Moment_Inertia*Omega) + U_torque);
Xdot(8:10) = v;
Xdot(11:13) = -mu_Earth/r_norm^3*r + Usrp; % + Udrag; % 

Xdot = Xdot(:);

end
% -------------------------------------------------------------------------
function u = SunTracker(P,k,I_c,r_target,v_target,r_chaser,v_chaser,q_chaser,Comega_chaser)

% Computing the reference frame
Nrhovec = r_target - r_chaser; n1 = [1 0 0]'; n3 = v_chaser/norm(v_chaser);
r3 = Nrhovec/norm(Nrhovec);

if norm(cross(r3,n3)) ~= 0
    r1 = cross(r3,n3)/norm(cross(r3,n3));
else
    if dot(r3,n3)<0
        r1 = cross(r3,n1)/norm(cross(r3,-n1));
    else
        r1 = cross(r3,n1)/norm(cross(r3,n1));
    end
end

r2 = cross(r3,r1)/norm(cross(r3,r1));
RN = [r1'; r2';r3'];

% Computing the reference angular velocity, angular velocity error,
% and the quaternion error
CN = Quaternion_to_DCM(q_chaser);
CR = CN*RN';
Nrhodotvec = v_target - v_chaser; rho = norm(Nrhovec);
omega_check = cross(Nrhovec,Nrhodotvec)/rho^2;
COmega_r = CN*omega_check;
CdelOmega = Comega_chaser - COmega_r; % Angular velocity error
q_errorr = DCM2Quaternion(CR); epsilon = q_errorr(2:4,:); % Quaternion error


% Computing the angular acceleration of the reference frame
COmega_tilde=[0 -Comega_chaser(3) Comega_chaser(2);...
    Comega_chaser(3) 0 -Comega_chaser(1);...
    -Comega_chaser(2) Comega_chaser(1) 0];
COmegadot_r = -I_c\(COmega_tilde*I_c*COmega_r); % Angular acceleratio of the reference frame

% Implementing the controller
u = - k*epsilon - P*CdelOmega ...
    + I_c*(COmegadot_r - COmega_tilde*COmega_r)...
    + COmega_tilde*I_c*Comega_chaser;
end
% -------------------------------------------------------------------------
function u = VelocityTracker(P,k,I_c,r_target,v_target,r_chaser,v_chaser,q_chaser,Comega_chaser)

% Computing the reference frame
T = v_chaser/norm(v_chaser);
W = cross(r_target,v_chaser)/norm(cross(r_target,v_chaser));
N = cross(T,W)/norm(cross(T,W));
M2 = [cos(pi) 0 -sin(pi);
      0       1        0;
      sin(pi) 0 cos(pi)];
Nrhovec = r_target - r_chaser; r3 = Nrhovec/norm(Nrhovec);
if dot(W,r3)>=0
    RN = [N'; T'; W'];
else
    RN = M2*[N'; T'; W'];
end
% Computing the reference angular velocity, angular velocity error,
% and the quaternion error
CN = Quaternion_to_DCM(q_chaser);
CR = CN*RN';
rho = norm(r_chaser);
omega_check = cross(r_chaser,v_chaser)/rho^2;
COmega_r = CN*omega_check;
CdelOmega = Comega_chaser - COmega_r; % Angular velocity error
q_errorr = DCM2Quaternion(CR); epsilon = q_errorr(2:4,:); % Quaternion error

% Computing the angular acceleration of the reference frame
COmega_tilde=[0 -Comega_chaser(3) Comega_chaser(2);...
    Comega_chaser(3) 0 -Comega_chaser(1);...
    -Comega_chaser(2) Comega_chaser(1) 0];
COmegadot_r = -I_c\(COmega_tilde*I_c*COmega_r); % Angular acceleratio of the reference frame

% Implementing the controller
u = - k*epsilon - P*CdelOmega ...
    + I_c*(COmegadot_r - COmega_tilde*COmega_r)...
    + COmega_tilde*I_c*Comega_chaser;
end

% -------------------------------------------------------------------------
function [a_drag] = Drag(r,v,q,Spacecraft)
format long
% Description of the sun location relative to the sailcraft body frame
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ECI_to_BodyFrame = Quaternion_to_DCM(q);
n1 = ECI_to_BodyFrame(1,:);
n2 = ECI_to_BodyFrame(2,:);
n3 = ECI_to_BodyFrame(3,:);
omega_Earth = [0 0 7.292124e-5]'; 
vrel = v - cross(omega_Earth,r);
vhat = vrel/norm(vrel);

s = Spacecraft.side; % Side of the cubesat
L = Spacecraft.L;    % length of the solarsail
w = Spacecraft.l;    % width of the solarsail
CD = Spacecraft.CD;    % width of the solarsail
mass = Spacecraft.M+Spacecraft.m; % Mass of the spacecraft (in Kg)

area_n1 = s^2*n1';
area_n2 = s^2*n2';
area_n3 = L*w*n3';
AreaMassRatio = (1/mass)*( abs(dot(area_n1,vhat)) ...
    + abs(dot(area_n2,vhat)) + abs(dot(area_n3,vhat))  );

rho = desity(r); % Earth's atmospheric desity
a_drag = -1/2*CD*AreaMassRatio*rho*dot(vrel,vrel)*vhat*1e3; % Acceleration on km/s^2

end
% -------------------------------------------------------------------------
function u = u_Cframe_QuaternionTracking(P,k,I_c,r_target,v_target,r_chaser,v_chaser,q_chaser,Comega_chaser)

% Computing the reference frame
Nrhovec = r_target - r_chaser; n1 = [1 0 0]'; n3 = [0 0 1]';
r3 = Nrhovec/norm(Nrhovec);
if norm(cross(r3,n3)) ~= 0
    r1 = cross(r3,n3)/norm(cross(r3,n3));
else
    if dot(r3,n3)<0
        r1 = cross(r3,n1)/norm(cross(r3,-n1));
    else
        r1 = cross(r3,n1)/norm(cross(r3,n1));
    end
end
r2 = cross(r3,r1)/norm(cross(r3,r1));
RN = [r1'; r2';r3'];

% Computing the reference angular velocity, angular velocity error,
% and the quaternion error
CN = Quaternion_to_DCM(q_chaser);
CR = CN*RN';
Nrhodotvec = v_target - v_chaser; rho = norm(Nrhovec);
omega_check = cross(Nrhovec,Nrhodotvec)/rho^2;
COmega_r = CN*omega_check;
CdelOmega = Comega_chaser - COmega_r; % Angular velocity error
q_errorr = DCM2Quaternion(CR); epsilon = q_errorr(2:4,:); % Quaternion error


% Computing the angular acceleration of the reference frame
COmega_tilde=[0 -Comega_chaser(3) Comega_chaser(2);...
    Comega_chaser(3) 0 -Comega_chaser(1);...
    -Comega_chaser(2) Comega_chaser(1) 0];
COmegadot_r = -I_c\(COmega_tilde*I_c*COmega_r); % Angular acceleratio of the reference frame

% Implementing the controller
u = - k*epsilon - P*CdelOmega ...
    + I_c*(COmegadot_r - COmega_tilde*COmega_r)...
    + COmega_tilde*I_c*Comega_chaser;
end
% -------------------------------------------------------------------------
function [a_srp,XS,VS] = FlatPlate_srp(time,X,Spacecraft,q)
format long
global JD

% Description of the sun location relative to the sailcraft body frame
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ECI_to_BodyFrame = Quaternion_to_DCM(q);

JulianDay = JD+time/86400;
[R_Earth, V_Earth, ~] = Ephem(JulianDay,3,'EME2000'); % Earth position and velocity measured from the Sun
XS = -R_Earth; % Sun's position measured from the Earth's ECI frame
VS = -V_Earth; % Sun's velocity measured from the Earth's ECI frame
Xrel = XS - X(1:3); NomrDelR = (Xrel.'*Xrel).^(1/2);
s = Xrel/NomrDelR; s = ECI_to_BodyFrame*s;

% Collecting sailcraft characteristics
rhos = Spacecraft.rhos;rhod = Spacecraft.rhod;
area = (Spacecraft.L*Spacecraft.l); % Area of the plate (in m^2)
Cr   = Spacecraft.Cr;
mass = Spacecraft.M+Spacecraft.m; % Mass of the spacecraft (in Kg)
Pcrp = 1357/(299792458); % Solar pressure (in﻿[N/m^2])
AU = 149597870.7; % Distance between Earth and Sun (in Km)
rs_squared = norm(Xrel/AU)^2;
P = Cr*Pcrp*area/(mass*rs_squared)*1e-3; % in km/s^2
n = [0 0 1]'; Costheta = n.'*s;

% Computing the external force exerted on the sailcraft (in km/s^2)
if isa(Costheta,'casadi.SX')
    a_srp  = if_else(Costheta<0, zeros(3,1),...
        -P*Costheta*((1-rhos)*s+(rhod+2*rhos*Costheta)*n));
else
    if Costheta<0
        a_srp = zeros(3,1);
    else
        a_srp = -P*Costheta*((1-rhos)*s+(rhod+2*rhos*Costheta)*n);
    end
end

end
% -------------------------------------------------------------------------
function [rho,T,v] = desity(X)
% Reference: https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
global ae
r = norm(X(1:3));
h = (r - ae)*1e3; % altitude above the earth surface (in meters)
if isa(h,'casadi.SX')
    temperature  = if_else(h<11000, 15.04-0.00649*h,...
        if_else(h>11000 & h<25000, -56.46,-131.21+0.00299*h));
    pressure = if_else(h<11000, 101.29*((temperature+273.1)/288.08)^5.256,...
        if_else(h>11000 & h<25000,22.65.*exp(1.73-0.000157*h),...
        2.488*((temperature+273.1)/216.6)^(-11.388)));
else
    if h<11000 % Troposphere
        temperature = 15.04-0.00649*h;
        pressure    = 101.29*((temperature+273.1)/288.08)^5.256;
    elseif  h>11000 && h<25000 % Law Stratosphre
        temperature = -56.46;
        pressure    = 22.65.*exp(1.73-0.000157*h);
    else % Upper Stratosphre
        temperature = -131.21+0.00299*h;
        pressure    = 2.488*((temperature+273.1)/216.6)^(-11.388);
    end
end
rho = pressure./(0.2869*(temperature+273.1)); % Atmospheric density
T   = temperature+273.1; % Temperature
v   = sqrt(1.4*286*T); % Speed of the sound

end
% -------------------------------------------------------------------------
function C = Quaternion_to_DCM(x)

if isa(x,'casadi.SX')
    C = [x(1)^2+x(2)^2-x(3)^2-x(4)^2 2*(x(2)*x(3)+x(1)*x(4)) 2*(x(2)*x(4)-x(1)*x(3));...
        2*(x(2)*x(3)-x(1)*x(4)) x(1)^2-x(2)^2+x(3)^2-x(4)^2 2*(x(3)*x(4)+x(1)*x(2));...
        2*(x(2)*x(4)+x(1)*x(3)) 2*(x(3)*x(4)-x(1)*x(2)) x(1)^2-x(2)^2-x(3)^2+x(4)^2];
    
else
    if (size(x,1) ~= 4)
        error('Quaternion parameter set must be 4-by-n matrix.')
    end
    if (abs(norm(x)-1)> 0.000001)
        error('Quaternions must have unit norm.')
    end
    C = [x(1)^2+x(2)^2-x(3)^2-x(4)^2 2*(x(2)*x(3)+x(1)*x(4)) 2*(x(2)*x(4)-x(1)*x(3));...
        2*(x(2)*x(3)-x(1)*x(4)) x(1)^2-x(2)^2+x(3)^2-x(4)^2 2*(x(3)*x(4)+x(1)*x(2));...
        2*(x(2)*x(4)+x(1)*x(3)) 2*(x(3)*x(4)-x(1)*x(2)) x(1)^2-x(2)^2-x(3)^2+x(4)^2];
end


end
% -------------------------------------------------------------------------
function b = DCM2Quaternion(C)
% Converts DCM to quaternions via the Stanley Method

B2(1) = (1+trace(C))/4;
B2(2) = (1+2*C(1,1)-trace(C))/4;
B2(3) = (1+2*C(2,2)-trace(C))/4;
B2(4) = (1+2*C(3,3)-trace(C))/4;

[~,i] = max(B2);
switch i
    case 1
        b(1) = sqrt(B2(1));
        b(2) = (C(2,3)-C(3,2))/4/b(1);
        b(3) = (C(3,1)-C(1,3))/4/b(1);
        b(4) = (C(1,2)-C(2,1))/4/b(1);
    case 2
        b(2) = sqrt(B2(2));
        b(1) = (C(2,3)-C(3,2))/4/b(2);
        if (b(1)<0)
            b(2) = -b(2);
            b(1) = -b(1);
        end
        b(3) = (C(1,2)+C(2,1))/4/b(2);
        b(4) = (C(3,1)+C(1,3))/4/b(2);
    case 3
        b(3) = sqrt(B2(3));
        b(1) = (C(3,1)-C(1,3))/4/b(3);
        if (b(1)<0)
            b(3) = -b(3);
            b(1) = -b(1);
        end
        b(2) = (C(1,2)+C(2,1))/4/b(3);
        b(4) = (C(2,3)+C(3,2))/4/b(3);
    case 4
        b(4) = sqrt(B2(4));
        b(1) = (C(1,2)-C(2,1))/4/b(4);
        if (b(1)<0)
            b(4) = -b(4);
            b(1) = -b(1);
        end
        b(2) = (C(3,1)+C(1,3))/4/b(4);
        b(3) = (C(2,3)+C(3,2))/4/b(4);
end
b = b';
end
% -------------------------------------------------------------------------
function Xdot = SRPNotCtrlODE(t, X, Target)
global mu_Earth
q = X(1:4)/norm(X(1:4));
Omega = X(5:7);
r = X(8:10); v = X(11:13);
r_norm = norm(r);
tilde = @(v) [ 0     -v(3)  v(2) ;
    v(3)   0    -v(1) ;
    -v(2)   v(1)  0   ];

% Rotational state differential equation
Moment_Inertia = diag(Target.Moment_Of_Inertia_Calculator());
W_Matrix       = [0, -Omega(1), -Omega(2), -Omega(3);
    Omega(1),  0,  Omega(3), -Omega(2);
    Omega(2), -Omega(3),  0,  Omega(1);
    Omega(3), Omega(2), -Omega(1), 0];

% Computing SRP acceleration (Flat Plate model)
[a_srp,~,~] = FlatPlate_srp(t, r, Target, q);
BN = Quaternion_to_DCM(q);
Usrp = BN\a_srp;

% Integrating the the target trajectory
Xdot(1:4) = 1/2*W_Matrix*q;
Xdot(5:7) = Moment_Inertia\(-tilde(Omega)*(Moment_Inertia*Omega));
Xdot(8:10) = v;
Xdot(11:13) = -mu_Earth/r_norm^3*r + Usrp;
Xdot = Xdot(:);
end
% -------------------------------------------------------------------------
function Xdot = SRPandDragODE(t, X, Target)
global mu_Earth
q = X(1:4)/norm(X(1:4));
Omega = X(5:7);
r = X(8:10); v = X(11:13);
q2 = X(14:17)/norm(X(14:17));
Omega2 = X(18:20);

r_norm = norm(r);
tilde = @(v) [ 0     -v(3)  v(2) ;
    v(3)   0    -v(1) ;
    -v(2)   v(1)  0   ];
BigTilde = @(v) [0, -v(1), -v(2), -v(3);
                v(1),  0,  v(3), -v(2);
                v(2), -v(3),  0,  v(1);
                v(3), v(2), -v(1), 0 ];


Moment_Inertia = diag(Target.Moment_Of_Inertia_Calculator());
% ======================================================================= %
% Switching condition for sun tracker
[u_srp, XS, VS] = FlatPlate_srp(t, r, Target, q);
adot = da_dt(r,v,q2,t);
% ======================================================================= %
P = 30*diag([1 2 3]); k = 30; % Controller gains
if adot<0 % Be align with the velocity vector
    U_torque = VelocityTracker(P,k,Moment_Inertia,XS,VS,r,v,q,Omega);
else
    U_torque = SunTracker(P,k,Moment_Inertia,XS,VS,r,v,q,Omega);
end

% Compting the the acceleration on the plate
BN = Quaternion_to_DCM(q);
Usrp = BN\u_srp;
Udrag = Drag(r,v,q,Target);
% DragNorm = norm(Udrag);
% Integrating the the target trajectory
Xdot(1:4) = 1/2*BigTilde(Omega)*q;
Xdot(5:7) = Moment_Inertia\(-tilde(Omega)*(Moment_Inertia*Omega) + U_torque);
Xdot(8:10) = v;
Xdot(11:13) = -mu_Earth/r_norm^3*r + Usrp + Udrag;

% ======================================================================= %
% Compting the torque required to track the sun
U_torque2   = SunTracker(P,k,Moment_Inertia,XS,VS,r,v,q2,Omega2);
Xdot(14:17) = 1/2*BigTilde(Omega2)*q2;
Xdot(18:20) = Moment_Inertia\(-tilde(Omega2)*(Moment_Inertia*Omega2) + U_torque2);


Xdot = Xdot(:);

end

