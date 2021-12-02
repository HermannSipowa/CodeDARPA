function [Xdot, DelUsrp, U_eci_Chasser, U_eci_Target,Costheta,del_ag] = ControlsAccelerationCheck(t,X,Target,Chasser,KMat,PMat,AgentNom)
format long
global mu_Earth Tc
ti = t*Tc;

tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];

X_chief = full([Chief_x(tt);  Chief_y(tt);  Chief_z(tt);...
            Chief_dx(tt); Chief_dy(tt); Chief_dz(tt)]);       

% Collecting relevant quantities    
q_Chasser = X(1:4)/norm(X(1:4)); OmegaBody = X(5:7);
r_target = X_chief(1:3); v_target = X_chief(4:6);
rt_norm = norm(r_target);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

rho = X(1:3);  rho_prime = X(4:6);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
r_chasser = TN\rc; % ECI position of the chasser

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
% F_body = F_CanonBall(ti,r_chasser,Chasser); % SRP force on the Chasser
[F_body,Costheta] = F_SPR_FlatPlate(ti, r_chasser, Chasser, q_Chasser);


%***************************************%
% Rotational state differential equation
I_Chasser  = diag(Chasser.Moment_Of_Inertia_Calculator());
Thau       = LyapunovAttitudeTrackingCtrl(t,X,KMat,PMat,AgentNom,I_Chasser);   
Omega_matrix_Chasser = [0, -OmegaBody(1), -OmegaBody(2), -OmegaBody(3);
                        OmegaBody(1),  0,  OmegaBody(3), -OmegaBody(2);
                        OmegaBody(2), -OmegaBody(3),  0,  OmegaBody(1);
                        OmegaBody(3), OmegaBody(2), -OmegaBody(1),  0];

qdot          = 1/2*Omega_matrix_Chasser*q_Chasser;
OmegaBody_dot = I_Chasser\(-tilde(OmegaBody)*(I_Chasser*OmegaBody)+Thau);

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
Coriolis_Effect = - 2*tilde(Omega)*rho_prime ...
                  - tilde(Omega)*(tilde(Omega)*rho) ...
                  - tilde(Omegadot)*rho;


% Integrating the relative trajectory
Rho_dot = [qdot; OmegaBody_dot; rho_prime; DelUsrp + del_ag + Coriolis_Effect];

% Collecting the rates of change
Xdot = Rho_dot*Tc;

end