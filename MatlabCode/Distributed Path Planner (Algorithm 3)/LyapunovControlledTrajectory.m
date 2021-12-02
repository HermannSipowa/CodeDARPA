function [Xdot] = LyapunovControlledTrajectory(t,X_Chasser,KMat,PMat,G,AgentID,AgentNom)
global Tc
u = LyapunovTrackingCtrl(t,X_Chasser,KMat,PMat,G,AgentID,AgentNom);
f_deputy  = RelativeDynamics(t,X_Chasser);
%% Integrating the deputy's trajectory
Xdot(1:3,1) = f_deputy(1:3);
Xdot(4:6,1) = f_deputy(4:6) + u;
Xdot = Xdot*Tc;
end