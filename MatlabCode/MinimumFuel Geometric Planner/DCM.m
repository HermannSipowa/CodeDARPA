function [CN] = DCM(r,v)
%--------------------------------------------------------------------------
% [CN] = DCM(Position1,Velocity1)
%
% Relative orientation of the chief relative to the inertial frame
%
% [inputs]: - [r]  : position vector (in Km)
%           - [v]  : Velocity vector (in Km/s)
%
% [outputs]: -[CN] : Direct cosine matrix from ECI to LVLH
%--------------------------------------------------------------------------
tilde = @(v) [0    -v(3)  v(2);
              v(3)   0   -v(1);
             -v(2)  v(1)    0];
         
Cross_rv = tilde(r)*v;
Norm_r = (r.'*r).^(1/2); 
NormCross_rv = (Cross_rv.'*Cross_rv).^(1/2);

x_hat = r/Norm_r;
z_hat = Cross_rv/NormCross_rv; 

Cross_zxhats = tilde(z_hat)*x_hat;
NormCross_zxhats = (Cross_zxhats.'*Cross_zxhats).^(1/2);
y_hat = Cross_zxhats/NormCross_zxhats;


CN = [x_hat'; y_hat'; z_hat'];

end