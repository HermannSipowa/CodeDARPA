function [xdot] = TH_ode(t,x)
format long
global e_chief
k = 1+e_chief*cos(t);
A = [zeros(3) eye(3);
    3/k  0  0   0   2  0;
    0    0  0   -2  0  0;
    0    0  -1  0   0  0];
U1 = [u11(t) u12(t) u13(t)]';
U2 = [u21(t) u22(t) u23(t)]';
u = [U1 U2]; 
B = [zeros(3); eye(3)];
[l,~] = size(x);
x = reshape(x,6,l/6);
xdot = A*x + B*u;
xdot = xdot(:);
end

