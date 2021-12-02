function [rho,T,v] = desity(X)
% Reference: https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
global ae
r = norm(X(1:3));
h = (r - ae)*1e3; % altitude above the earth surface (in meters)
if isa(h,'casadi.SX')
    temperature  = if_else(h<11000, 15.04-0.00649*h,...
        if_else(h>11000 & h<25000, -56.46,-131.21+0.00299*h));
    pressure = if_else(h<11000, 101.29*((temperature+273.1)/288.08)^5.256,...
        if_else(h>11000 & h<25000,22.65.*exp(1.73-0.000157*h),...
        2.488*((temperature+273.1)/216.6)^(-11.388)));
else
    if h<11000 % Troposphere
        temperature = 15.04-0.00649*h;
        pressure    = 101.29*((temperature+273.1)/288.08)^5.256;
    elseif  h>11000 && h<25000 % Law Stratosphre
        temperature = -56.46;
        pressure    = 22.65.*exp(1.73-0.000157*h);
    else % Upper Stratosphre
        temperature = -131.21+0.00299*h;
        pressure    = 2.488*((temperature+273.1)/216.6)^(-11.388);
    end
end
rho = pressure./(0.2869*(temperature+273.1)); % Atmospheric density
T   = temperature+273.1; % Temperature
v   = sqrt(1.4*286*T); % Speed of the sound

end