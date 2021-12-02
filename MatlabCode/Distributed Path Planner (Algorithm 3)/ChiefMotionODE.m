function Xdot = ChiefMotionODE(t, X, Target)
format long
global mu_Earth

r_target = X(1:3,:); v_target = X(4:6,:);
rt_norm = norm(r_target);

%% Computing SRP acceleration (Cannonball model)
[U_eci_Target] = F_CanonBall(t,r_target,Target); % SRP force on the Target

%% Integrating the the target trajectory
Xdot(1:3,:) = v_target;
Xdot(4:6,:) = -mu_Earth/rt_norm^3*r_target + U_eci_Target;

end