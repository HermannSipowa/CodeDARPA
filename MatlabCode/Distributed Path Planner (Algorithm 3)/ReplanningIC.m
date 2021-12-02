% --------------------- Initial condition for pdepe ----------------------
function u0 = ReplanningIC(x,AgentNom,AgentID) % Initial condition for AGHF  X0, Xf, To, T,

u0 = AgentNom(AgentID).DesiredStateInterpolant(x)';

end