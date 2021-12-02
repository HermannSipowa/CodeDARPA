function [Xdot] = NL_ode(t,x,agent)
global mu_Earth Rrho Tc

format long
X_chief = full([Chief_x(t);  Chief_y(t);  Chief_z(t);...
    Chief_dx(t); Chief_dy(t); Chief_dz(t)]);
if agent == 1
    u = [u11(t) u12(t) u13(t)]';
elseif agent == 2
    u = [u21(t) u22(t) u23(t)]';
end

tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

% Collecting relevant quantities
rt = X_chief(1:3); vt = X_chief(4:6); rt_norm = (rt.'*rt).^(1/2);
rho_bar = x(1:3,:);  drho_bar = x(4:6,:);

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

B = [zeros(3); eye(3)];
Xdot = [drho_bar; (del_ag + Coriolis_Effect)];
Xdot = (blkdiag(Rrho,Rrho)\Xdot)*Tc + B*u;

end