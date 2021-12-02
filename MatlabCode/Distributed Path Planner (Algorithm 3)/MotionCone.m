function [s] = MotionCone(rho,AnglePhi,R,Origin,Color)
hold on
alpa = 0.1;
if isempty(rho)
    rho=2;
end
if isempty(AnglePhi)
    AnglePhi=pi/6;
end
if isempty(R)
    R=eye(3);
end
if isempty(Origin)
    Origin = [0;0;0];
end
if isempty(Color)
    Color = [145,90,7]/255;
end



AngleGamma = pi/2-AnglePhi;
AngleTheta = 2*pi;
theta=linspace(0,AngleTheta,20);
phi=linspace(0,AnglePhi,20);
[THETA,PHI]=meshgrid(theta,phi);

X=rho*sin(PHI).*cos(THETA);
Y=rho*sin(PHI).*sin(THETA);
Z=rho*cos(PHI);

XYZ = [Y(:) X(:) Z(:)];
RotXYZ = XYZ*R';
Xqr = Origin(1)+reshape(RotXYZ(:,1),size(X,1), []);
Yqr = Origin(2)+reshape(RotXYZ(:,2),size(Y,1), []);
Zqr = Origin(3)+reshape(RotXYZ(:,3),size(Z,1), []);
s = surf(Xqr,Yqr,Zqr);
s.FaceAlpha = alpa;
s.EdgeColor = 'none';
s.FaceColor = 'b';



gamma=linspace(0,AngleGamma,20);
[THETA,GAMMA]=meshgrid(theta,gamma);
X=rho*sin(PHI).*cos(THETA);
Y=rho*sin(PHI).*sin(THETA);
Z=rho.*sin(GAMMA);


XYZ = [X(:) Y(:) Z(:)];
RotXYZ = XYZ*R';
Xqr = Origin(1)+reshape(RotXYZ(:,1),size(X,1), []);
Yqr = Origin(2)+reshape(RotXYZ(:,2),size(Y,1), []);
Zqr = Origin(3)+reshape(RotXYZ(:,3),size(Z,1), []);
s = surf(Xqr,Yqr,Zqr);
s.FaceAlpha = alpa;
s.EdgeColor = 'none';
s.FaceColor = 'b';


xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
% axis equal

end