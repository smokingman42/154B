function [wx,wy,M,d,l,z] = get_wxwy(n,rho,V,AoA,Cd,CM,nz)
% to calculate wx wy under equilibrium 
b = 10.82;                               % m          full span
c = 1.5;                                 % m          chord length
S = b*c;                                 % m^2        surface area
AR = b^2/S;                              %            aspect ratio
e = 0.79;                                %            Oswald efficiency     
W = 1100*9.81;                                % N          weight
   
L = n*W; % N          lift
CL = 2*L/rho/V^2/S;                      %            lift coefficient
D = 0.5*rho*V^2*S*(Cd + CL^2/pi/AR/e);   % N          drag 

M = CM*.5*rho*V^2*S*c;
z = 0:b/2/nz:b/2;                        % root to tip (half span)

% lift distribution
l_rect = L/b.*ones(1,nz+1);
l_ellip = (4*L/pi/b).*sqrt(1-(2.*z./b).^2);
l = (l_rect + l_ellip)./2;               % N/m

% figure
% plot(z,l_ellip,z,l_rect,z,l,'Linewidth',2)
% xlabel('z (m)')
% ylabel('Lift Distribution (N/m)')
% legend('lift elliptic distribution','lift rectangular distribution','combined lift distribution','location','west')

% drag distribution
d = 1.1*D/b*ones(1,nz+1);
nz2 = round(nz*.8,0)-1;
d(1:nz2)=D/b;

% figure
% plot(z,d,'Linewidth',2)
% xlabel('z (m)')
% ylabel('Drag Distribution (N/m)')
% legend('Drag Distribution','location','west')


% rotate into x-y coordinate
wy = cos(AoA).*l + sin(AoA).*d;
wx = -sin(AoA).*l + cos(AoA).*d;

% Note: wx and wy are defined from root to tip

% figure
% plot(z,wy,z,wx,'linewidth',2)
% xlabel('z (m)')
% ylabel('Distributed Load (N/m)')
% legend('w_y','w_x')