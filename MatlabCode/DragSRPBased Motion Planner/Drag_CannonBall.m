function acc = Drag_CannonBall(r,v,Spacecraft)

omega_Earth = [0 0 7.292124e-5]'; 
vrel = v - cross(omega_Earth,r);
vhat = vrel/norm(vrel);
CD = Spacecraft.CD;    % width of the solarsail
mass = Spacecraft.M+Spacecraft.m; % Mass of the spacecraft (in Kg)
Area = (2*pi*Spacecraft.r^2); % Area of the plate (in m^2)
AreaMassRatio = Area/mass;
rho = desity(r); % Earth's atmospheric desity

acc = -1/2*CD*AreaMassRatio*rho*dot(vrel,vrel)*vhat*1e3; % Acceleration on km/s^2
end