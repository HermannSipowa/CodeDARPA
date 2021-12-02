function [Xdot] = DragSRP_FlatPlateControlledOde(t,X,Target,Chasser,KMat,PMat,AgentNom)
format long
global mu_Earth Tc
ti = t*Tc;
tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];

% Collecting relevant quantities
Xchief = full(ChiefState(t));
r_target = Xchief(1:3); v_target = Xchief(4:6);
rt_norm = (r_target.'*r_target).^(1/2);
TN = DCM(r_target,v_target); % DCM's between target's RIC frame and ECI frame

q_Chasser = X(1:4)/((X(1:4).'*X(1:4)).^(1/2));
Omega_Chasser_body = X(5:7);
rho = X(8:10); rho_prime = X(11:13);
rc = [rt_norm+rho(1) rho(2) rho(3)]'; % RIC position of the chasser
rc_norm = (rc.'*rc).^(1/2); 
r_Chasser = TN\rc; % ECI position of the chasser



%***************************************%
% Rotational state differential equation
I_Chasser        = diag(Chasser.Moment_Of_Inertia_Calculator());
OmegaMat_Chasser = [0, -Omega_Chasser_body(1), -Omega_Chasser_body(2), -Omega_Chasser_body(3);
                    Omega_Chasser_body(1),  0,  Omega_Chasser_body(3), -Omega_Chasser_body(2);
                    Omega_Chasser_body(2), -Omega_Chasser_body(3),  0,  Omega_Chasser_body(1);
                    Omega_Chasser_body(3), Omega_Chasser_body(2), -Omega_Chasser_body(1),  0];
Thau = LyapunovAttitudeTrackingCtrl(t,X,KMat,PMat,AgentNom,I_Chasser); 
qdot      = 1/2*OmegaMat_Chasser*q_Chasser;
Omega_dot = I_Chasser\(-tilde(Omega_Chasser_body)*(I_Chasser*Omega_Chasser_body) + Thau);

%***************************************%

% Computingthe Drag and SRP accelerations (Cannonball model)
[U_eci_Target, Xrel, Vrel, beta_chief] = F_CanonBall(ti,Xchief,Target); % SRP force on the Target
FDeputy_body = F_SPR_FlatPlate(ti, r_Chasser, Chasser, q_Chasser);


% Translational state differential equation
BN = Quaternion_to_DCM(q_Chasser);
U_eci_Chasser = BN\FDeputy_body(1:3);

% Computing the angular velocity and angular acceleration of the target frame
h_vec   = tilde(r_target)*v_target; h_norm = (h_vec.'*h_vec).^(1/2); u_acc = U_eci_Target;
eh      = h_vec/h_norm; h_dot = tilde(r_target)*u_acc; Xrel_norm = (Xrel.'*Xrel).^(1/2);
eh_dot  = h_dot/h_norm - dot(h_vec,h_dot)*h_vec/h_norm^3;
acc_dot = beta_chief*(Vrel/Xrel_norm^3 - 3*dot(Xrel,Vrel)*Xrel/Xrel_norm^5);

    
Omega    = TN*( h_vec/rt_norm^2 + dot(u_acc,eh)*r_target/h_norm );  
Omegadot = TN*( h_dot/rt_norm^2 ...
                - 2*dot(r_target,v_target)*h_vec/rt_norm^4 ...
                - dot(h_vec,h_dot)*dot(u_acc,eh)*r_target/h_norm^3 ...
                + (dot(acc_dot,eh) + dot(u_acc,eh_dot))*r_target/h_norm ...
                + dot(u_acc,eh)*v_target/h_norm );

% Differential Drag acceleration
Udrag_Target = Drag_CannonBall(r_target,v_target,Target);         
v_Chasser = v_target + TN\(rho_prime + cross(Omega,rho));
UDrag_Chasser = Drag(r_Chasser,v_Chasser,q_Chasser,Chasser);

% Computing the relative acceleration
Del_ag  =  -mu_Earth/rc_norm^3*rc + mu_Earth/rt_norm^2*[1 0 0]'; % (2-body) gravitatinal
DelUsrp = TN*(U_eci_Chasser - U_eci_Target); % Delta-SRP
DelDrag = TN*(UDrag_Chasser-Udrag_Target);
Rotation_Effect = - 2*tilde(Omega)*rho_prime ...
                  - tilde(Omega)*(tilde(Omega)*rho) ...
                  - tilde(Omegadot)*rho;


% Integrating the relative trajectory
Rho_dot = [rho_prime; DelDrag + DelUsrp + Del_ag + Rotation_Effect];

% Collecting the rates of change
Xdot = [qdot; Omega_dot; Rho_dot]*Tc;

end