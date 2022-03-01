function [a_srp,Costheta] = F_SPR_FlatPlate(time,X,Spacecraft,q)
% -------------------------------------------------------------------------
% [F_body] = F_SPR_FlatPlate(t, X, Spacecraft, Aug_X)
% 
% The SRP force model comes from the following paper:
% J. Van Der Ha and V. Lappas, "Long-Term Attitude Drift of Spinning 
% Spacecraft Under Solar Radiation Torques", Journal of Guidance, Control,
% and Dynamics, vol. 30, no. 5, pp. 1470-1479, 2007.
%
% [inputs]:  - [time]       : simulation time (for computing the sun sun ephemeris)
%            - [X]          : Spacecraft inertial state
%            - [Spacecraft] : class defining a solar sail
%
% [outputs]: - [F_body] : F(1:3,1) = Acceleration on acting on the flat plate,
%                         F(4:6,1) = Torque acting on the flat plate 
%                         expressed in the body fixed frame of the sailcraft
% -------------------------------------------------------------------------
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
a_srp = -P*Costheta*((1-rhos)*s+(rhod+2*rhos*Costheta)*n);
end

