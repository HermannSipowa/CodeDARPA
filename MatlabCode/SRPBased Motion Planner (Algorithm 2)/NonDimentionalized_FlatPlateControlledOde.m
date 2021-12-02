function [Xdot] = NonDimentionalized_FlatPlateControlledOde(t,X,Target,Chasser,KMat,PMat,AgentNom)
format long
global mu_Earth Tc
ti = t*Tc;

tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];
IMat = diag(Chasser.Moment_Of_Inertia_Calculator());
% Collecting relevant quantities
X_chief = full([Chief_x(t);  Chief_y(t);  Chief_z(t);...
                Chief_dx(t); Chief_dy(t); Chief_dz(t)]);
r_target = X_chief(1:3); v_target = X_chief(4:6);
rt_norm = norm(r_target);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

qBody = X(1:4)/norm(X(1:4));
OmegaBody = X(5:7);
rho = X(8:10); rho_prime = X(11:13);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
rChasser = TN\rc; % ECI position of the chasser

% Computing SRP acceleration (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,X_chief,Target); % SRP force on the Target
FDeputy_body = F_SPR_FlatPlate(ti, rChasser, Chasser, qBody);


%***************************************%
% Computing the attitude tracking control

% % Attitude and angular velocity of the reference frame
% XRef = AgentNom.DesiredStateInterpolant(t)'; 
% qRef = XRef(1:4); OmegaRef = XRef(5:7); 
% 
% % Angular acceleration of the reference frame
% XdotRef = AgentNom.FDesiredInterpolant(t)'; 
% OmegadotRef = XdotRef(5:7);
% qRef = qRef/norm(qRef);
% 
% % rotation matrix between body frame and the reference frame 
% BR = Quaternion_to_DCM(qBody)*Quaternion_to_DCM(qRef)'; 
% epsilon = DCM2Quaternion(BR);
% delOmega = OmegaBody-OmegaRef;
% 
% % Defining the attitude tracking control (for quaternions)
% Thau = - KMat*epsilon(2:end) - PMat*delOmega ...
%     + IMat*(OmegadotRef - tilde(OmegaBody)*OmegaRef)...
%     + tilde(OmegaBody)*IMat*OmegaBody;
% 
% check = u-Thau
%***************************************%
% Rotational state differential equation
OmegaBodyMat = [0, -OmegaBody(1), -OmegaBody(2), -OmegaBody(3);
                OmegaBody(1),  0,  OmegaBody(3), -OmegaBody(2);
                OmegaBody(2), -OmegaBody(3),  0,  OmegaBody(1);
                OmegaBody(3), OmegaBody(2), -OmegaBody(1),  0];
Thau = LyapunovAttitudeTrackingCtrl(t,X,KMat,PMat,AgentNom,IMat); 
qBodyDot     = 1/2*OmegaBodyMat*qBody;
OmegaBodyDot = IMat\(-tilde(OmegaBody)*(IMat*OmegaBody) + Thau);
%***************************************%
% Translational state differential equation
BN = Quaternion_to_DCM(qBody);
U_eci_Chasser = BN\FDeputy_body(1:3);

% Computing the angular velocity and angular acceleration of the target frame
h_vec   = tilde(r_target)*v_target; h_norm = norm(h_vec); u_acc = U_eci_Target;
eh      = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; Xrel_norm = norm(Xrel);
eh_dot  = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
acc_dot = beta_chief*(Vrel/Xrel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/Xrel_norm^5);

    
OmegaLVLH     = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
OmegadotLVLH  = TN*( h_dot/rt_norm^2 ...
                - 2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% Computing the relative acceleration
del_ag  =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelUsrp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP
Rotation_Effect = - 2*tilde(OmegaLVLH)*rho_prime ...
                  - tilde(OmegaLVLH)*(tilde(OmegaLVLH)*rho) ...
                  - tilde(OmegadotLVLH)*rho;


% Integrating the relative trajectory
Rho_dot = [rho_prime; DelUsrp + del_ag + Rotation_Effect];

% Collecting the rates of change
Xdot = [qBodyDot; OmegaBodyDot; Rho_dot]*Tc;

end