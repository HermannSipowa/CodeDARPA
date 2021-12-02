% --------------------------- PDE for pdepe ------------------------------
function [c,f,s] = mypdexpde(x,t,u,DuDx,k,N) % Define PDE; right-hand-side of AGHF
XChief = full(ChiefState(x));
LL = full(L(u,DuDx,k,x,XChief));
CasADiresult = full(dLdx(u,DuDx,k,x,XChief,LL))';
NN = length(u);
CasADir = [CasADiresult(1:NN), CasADiresult(NN+1:2*NN)];

f =  CasADir(:,2); 
s = -CasADir(:,1); 
c = ones(N,1);

end