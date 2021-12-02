function [u] = LyapunovAttitudeTrackingCtrl(t,X_Chasser,KMat,PMat,AgentNom,IMat)

tilde = @(v) [ 0     -v(3)  v(2) ;
               v(3)   0    -v(1) ;
              -v(2)   v(1)  0   ];
% Attitude and angular velocity of the body frame
q = X_Chasser(1:4); Omega = X_Chasser(5:7);

% Attitude and angular velocity of the reference frame
X_ref = AgentNom.DesiredStateInterpolant(t)'; 
q_ref = X_ref(1:4); Omega_r = X_ref(5:7); 

% Angular acceleration of the reference frame
Xdot_ref = AgentNom.FDesiredInterpolant(t)'; 
Omegadot_r = Xdot_ref(5:7);
q = q/norm(q); q_ref = q_ref/norm(q_ref);

% rotation matrix between body frame and the reference frame 
BN = Quaternion_to_DCM(q);
NR = Quaternion_to_DCM(q_ref)';
BR = BN*NR; 
epsilon = DCM2Quaternion(BR);
delOmega = Omega-Omega_r;

% Defining the attitude tracking control (for quaternions)
u = - KMat*epsilon(2:end) - PMat*delOmega ...
    + IMat*(Omegadot_r - tilde(Omega)*Omega_r)...
    + tilde(Omega)*IMat*Omega;
end
