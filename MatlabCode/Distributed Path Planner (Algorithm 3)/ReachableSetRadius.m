function DelRMat = ReachableSetRadius(t,dt,AgentNom,NumAgent)
% This function compute the radius of the reachable set, whxich is
% considered as a sphere
% t: time at the center of the sphere X(t)
% dt: time to go before reaching the edge of the sphere x(t+dt)
r = nan(6,NumAgent); dr = nan(6,NumAgent);
for j = 1:NumAgent
    r(:,j)  = AgentNom(j).StatesExtendedInterpolant(t)';
    dr(:,j) = AgentNom(j).StatesExtendedInterpolant(t+dt)';
end
DelRadius = dr - r;
DelRMat = vecnorm(DelRadius(1:3,:),2,1);

end