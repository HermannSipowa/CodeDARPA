function u0 = mypdexic(x,X0,Xf,T) % Initial condition for AGHF
% Sinusoidal initial condition
freq1 = pi/2 + 2*pi;
freq2 = pi + 2*pi;
u0 = X0*cos(freq1*x/T) ...
     + Xf*sin(freq1*x/T) ...
     - (X0+Xf)*sin(freq2*x/T);
u0(1:4) = u0(1:4)/norm(u0(1:4)); % imposing the quaternion constraint
end