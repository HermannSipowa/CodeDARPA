% ---------------------= PDE for Replanning Sequence -----===--------------
function [c,f,s] = ReplanningPDE(x,t,u,DuDx,k,N,ObstaclesList,kb,pb,AgentNom)

Nobs      = length(ObstaclesList);
XChief    = full(ChiefState(x));
Obstacles = zeros(3,9);
for j = 1:Nobs
    ObsID = ObstaclesList(j);
    ObstacleStates = AgentNom(ObsID).ObstaclePredictionInterpolant(x)';
    Obstacles(:,j) = ObstacleStates(1:3,:);
end

[Lr1,Lr2,Lr3,Lr4,Lr5,Lr6,Lr7,Lr8,Lr9] = Lreplan(u,DuDx,k,x,XChief,Obstacles,kb,pb);
DL = Lreplandx(u,DuDx,k,x,XChief,Obstacles,kb,pb,Lr1,Lr2,Lr3,Lr4,Lr5,Lr6,Lr7,Lr8,Lr9);

if Nobs == 1
    CasADir = [full(DL(1,1:6))', full(DL(1,7:12))'];
elseif Nobs == 2
    CasADir = [full(DL(2,1:6))', full(DL(2,7:12))'];
elseif Nobs == 3
    CasADir = [full(DL(3,1:6))', full(DL(3,7:12))'];
elseif Nobs == 4
    CasADir = [full(DL(4,1:6))', full(DL(4,7:12))'];
elseif Nobs == 5
    CasADir = [full(DL(5,1:6))', full(DL(5,7:12))'];
elseif Nobs == 6
    CasADir = [full(DL(6,1:6))', full(DL(6,7:12))'];
elseif Nobs == 7
    CasADir = [full(DL(7,1:6))', full(DL(7,7:12))'];
elseif Nobs == 8
    CasADir = [full(DL(8,1:6))', full(DL(8,7:12))'];
elseif Nobs == 9
    CasADir = [full(DL(9,1:6))', full(DL(9,7:12))'];
end
    
% dL/dx
pLx  = CasADir(:,1); % EL(:,1); %
% dL/d(dot_x)
pLxd = CasADir(:,2); % EL(:,2); %

f = pLxd;
s = -pLx;
c = ones(N,1);

end % Define PDE; right-hand-side of AGHF