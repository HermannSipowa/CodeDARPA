function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t,X0,Xf,N) % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end