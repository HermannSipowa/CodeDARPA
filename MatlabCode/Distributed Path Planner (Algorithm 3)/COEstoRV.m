function [Position,Velocity,f] = COEstoRV(COE,mu)
%--------------------------------------------------------------------------
% [Position,Velocity] = COEstoRV(OE,mu)
%
% This function converts classical orbital elements to ECI cartesian
% coordinates. All angles here are in radians
%
% [inputs]:  - [COE] : Keplerian (osculating) orbital elements
%            - [mu]  : Body's gravitational constant
%
% [outputs]: - [Position]  : ECI position vector (in Km)
%            - [Velocity]  : ECI Velocity vector (in Km/s)
%--------------------------------------------------------------------------

a = COE(1); e = COE(2); inc = COE(3); W = COE(4); w = COE(5); M = COE(6);

%% Solving Kepler Equation
E = M; % Initial guess for the Eccentric Anomaly
if e < 1 % Elliptical case
    M_estimated = E - e*sin(E); % Calculated Mean motion
    Dm_dE = 1 - e*cos(E);
elseif e > 1 % Hyperbolic case
    M_estimated = e*sinh(E) - E;
    Dm_dE = e*cosh(E) - 1;
end

error = abs(M - M_estimated);
tol = 10^-15;

while error>tol
E = E + (M - M_estimated)/Dm_dE; % Eccentric anomaly

if e < 1 % Elliptical case
    M_estimated = E - e*sin(E);
    Dm_dE = 1 - e*cos(E);
elseif e > 1 % Hyperbolic case
    M_estimated = e*sinh(E) - E;
    Dm_dE = e*cosh(E) - 1;
end

error = abs(M - M_estimated);
end

%% Solving for the true anomaly 
if e < 1 % Elliptical case
    f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
elseif e > 1 % Hyperbolic case
    f = 2*atan(sqrt((1+e)/(e-1))*tanh(E/2));
end

%% Computing the magnetude of the position r
r = a*(1-e^2)/(1+e*cos(f));
rvec = [r*cos(f); r*sin(f); 0];
vvec = sqrt(mu/(a*(1-e^2)))*[-sin(f); e+cos(f); 0];
    
%% Computing the rotation matrices
R1=[cos(-W) sin(-W) 0; -sin(-W) cos(-W) 0; 0 0 1];
R2=[1 0 0; 0 cos(-inc) sin(-inc); 0 -sin(-inc) cos(-inc)];
R3=[cos(-w) sin(-w) 0; -sin(-w) cos(-w) 0; 0 0 1];

%% Computing the reverse rotation Matrix (i.e from the perifocal frame to the inertial frame)
DCM=R1*R2*R3;

%% Computing the vectors position and velocity in the inertial frame
X=DCM*rvec;
V=DCM*vvec;

%% Computing the position and velocity in the inertial frame
Position=[X(1,1) X(2,1)  X(3,1)]';
Velocity=[V(1,1) V(2,1) V(3,1)]';

end