function [Xdot] = RelativeDynamics(t,X)
global mu_Earth 
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
Xdot = XdotMat(:);

end