function Traj=Trajectory_GoldenAngle_GROG(ntviews,nx);
%Generate 2D golden-angle radial trajectory

% a=111.246117975;
Gn = (1 + sqrt(5))/2;
a=180/Gn;

radian=mod((0:a:(ntviews-1)*a)*pi/180,2*pi);
Rho=[-floor(nx/2):floor(nx/2)];
Rho=Rho(1:nx)+0.5;
% Rho=linspace(-(nx/2),nx/2,nx);

for jj=1:size(radian,2)
    X(:,jj)=Rho*sin(radian(1,jj));
    Y(:,jj)=-Rho*cos(radian(1,jj));
end

Traj=X+i*Y;
return