function [u,LuapCost,AGHFCost] = LyapunovTrackingCtrl(t,X_Chasser,KMat,PMat,G,AgentID,AgentNom)
global Tc

X_Desired = AgentNom(AgentID).DesiredStateInterpolant(t)';
f_desired = (1/Tc)*AgentNom(AgentID).FDesiredInterpolant(t)';

%% Computing the globally asymptotically stabilizing controller
f_deputy  = RelativeDynamics(t,X_Chasser);
del_r = X_Chasser(1:3) - X_Desired(1:3);
del_v = X_Chasser(4:6) - X_Desired(4:6);
u = f_desired(4:6) - f_deputy(4:6) - KMat*del_r - PMat*del_v;
xdot = [f_deputy(1:3); f_deputy(4:6)+u];
Error = xdot-f_desired;

%% Computing the cost function(s) to be minimize
LuapCost = 1/2*dot(del_v,del_v) + 1/2*del_r.'*KMat*del_r; % Lyapunov Cost function to minimize
AGHFCost = Error.'*G*Error; % Riemannian Cost function to minimize

end

