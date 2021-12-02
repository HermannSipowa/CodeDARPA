function [a_drag] = Drag(r,v,q,Spacecraft)
format long
% Description of the sun location relative to the sailcraft body frame
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ECI_to_BodyFrame = Quaternion_to_DCM(q);
n1 = ECI_to_BodyFrame(1,:);
n2 = ECI_to_BodyFrame(2,:);
n3 = ECI_to_BodyFrame(3,:);
omega_Earth = [0 0 7.292124e-5]'; 
vrel = v - cross(omega_Earth,r);
vhat = vrel/((vrel'*vrel)^(1/2));

s    = Spacecraft.side; % Side of the cubesat
L    = Spacecraft.L;    % length of the solarsail
w    = Spacecraft.l;    % width of the solarsail
CD   = Spacecraft.CD;  % width of the solarsail
mass = Spacecraft.M+Spacecraft.m; % Mass of the spacecraft (in Kg)

area_n1 = s^2*n1';
area_n2 = s^2*n2';
area_n3 = L*w*n3';
AreaMassRatio = (1/mass)*( abs(area_n1'*vhat) ...
    + abs(area_n2'*vhat) + abs(area_n3'*vhat)  );

rho = desity(r); % Earth's atmospheric desity
a_drag = -1/2*CD*AreaMassRatio*rho*(vrel'*vrel)*vhat*1e3; % Acceleration on km/s^2

end