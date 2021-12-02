function Xdot = DragSRPChiefMotionODE(t, X, Target)
format long
global mu_Earth

r = X(1:3,:); v = X(4:6,:); 

%% Computing SRP acceleration (Cannonball model)
USRP     = F_CanonBall(t,X,Target); % SRP force on the Target
Udrag    = Drag_CannonBall(r,v,Target);
Ugravity = -mu_Earth/(norm(r)^3)*r;
%% Integrating the the target trajectory
Xdot(1:3,:) = v;
Xdot(4:6,:) = Ugravity + USRP + Udrag;

end