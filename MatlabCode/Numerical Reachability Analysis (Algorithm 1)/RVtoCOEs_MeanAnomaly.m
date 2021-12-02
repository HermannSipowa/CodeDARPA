function [OBE] = RVtoCOEs_MeanAnomaly(r,v,nu)
%--------------------------------------------------------------------------
% [OBE]=RVtoCOEs(r,v,nu)
%
% This function coverts cartesian corrdinates into classical orbital
% elements (in radian)
% [inputs]: - [r]  : position vector (in Km)
%           - [v]  : Velocity vector (in Km/s)
%           - [nu] : Body's gravitational constant
%
% [outputs]: -[OBE] : Classical Obital Elements (angles in radian)
%--------------------------------------------------------------------------

k = [0 0 1]';    % normal vector in the k-direction
h = cross(r,v); % Anglar Momentum Vector
n = cross(k,h); % Node Vector
e_vec = (1/nu)*((norm(v)^2-nu/norm(r))*r-dot(r,v)*v); % Eccentricity Vector


P = norm(h)^2/nu; % Semi-Latus Rectum

e = norm(e_vec); % Eccentricity

i = acos(h(3)/norm(h)); % Inclination

Big_Omega = acos(n(1)/norm(n)); % Longitude of ascending node

if isa(Big_Omega,'casadi.SX')
    Big_Omega = if_else(n(2)<0,2*pi-Big_Omega,Big_Omega);
    
else
    if n(2)<0
        Big_Omega = 2*pi-Big_Omega;
    end
end

Little_Omega = acos(dot(n,e_vec)/(norm(n)*norm(e_vec))); % Argument of periapsis
if isa(Little_Omega,'casadi.SX')
    Little_Omega = if_else(e_vec(3)<0,2*pi-Little_Omega,Little_Omega);
else
    if e_vec(3)<0
        Little_Omega = 2*pi-Little_Omega;
    end
end


Nu = acos(dot(e_vec,r)/(norm(e_vec)*norm(r))); % True Anomaly
if isa(Nu,'casadi.SX')
    Nu = if_else(dot(r,v)<0,2*pi-Nu,Nu);
else
    Nu = real(Nu);
    if dot(r,v)<0
        Nu = 2*pi-Nu;
    end
end

E = 2*atan(sqrt((1-e)/(1+e))*tan(Nu/2));
M = E - e*sin(E);

a = P/(1-norm(e_vec)^2); % Semi-major Axis

OBE = [a; e; i; Big_Omega; Little_Omega; M];
end