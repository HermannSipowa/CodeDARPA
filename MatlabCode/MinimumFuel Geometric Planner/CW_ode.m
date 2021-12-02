function [xdot] = CW_ode(t,x)
% The EOM are nondimensionalization!!
format long 
global mu_dot Tc

% Computing the CW control and dynamic matrices
B = [zeros(3); eye(3)];
A  = [zeros(3) eye(3);
    3*mu_dot^2 0 0 0 2*mu_dot 0;
    0 0 0 -2*mu_dot 0 0;
    0 0 -mu_dot^2 0 0 0]*Tc;
U1 = [u11(t) u12(t) u13(t)]';
U2 = [u21(t) u22(t) u23(t)]';
u = [U1 U2]; 
% Reshaping the input into an 6 by n matrix
[l,~] = size(x);
x = reshape(x,6,l/6);
xdot = A*x + B*u;
xdot = xdot(:);
end