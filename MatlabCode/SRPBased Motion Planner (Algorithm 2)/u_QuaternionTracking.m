function u = u_QuaternionTracking(P,k,Omega,Omega_r,Omegadot_r,delOmega,epsilon,I)
tilde = @(v) [0    -v(3)  v(2);
    v(3)   0   -v(1);
    -v(2)  v(1)    0];

u = - k*epsilon -P*delOmega ...
    + I*(Omegadot_r - tilde(Omega)*Omega_r)...
    + tilde(Omega)*I*Omega;

end