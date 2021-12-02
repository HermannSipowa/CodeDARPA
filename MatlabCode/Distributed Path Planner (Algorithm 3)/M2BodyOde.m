function [Xdot] = M2BodyOde(t,X,mu_Earth)

Xdot(1:3,1)  = X(4:6,1);
Xdot(4:6,1) = -mu_Earth/norm(X(1:3,:))^3 * X(1:3,1);

end