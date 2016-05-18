clc; clear; close all;

%structure dimensions
b = 10.82;
c = 1.5; %m
kt = 0.001016; %m skin thickness
sh = 0.08; %m
st = 0.0025; %m spar thickness
bl = 0.015; %m bracket height
bt = 0.0025; %m bracket thickness`
theta = 30; %degrees
alpha = 10; %degrees
E = 70e9;  
G = E/(2+.6666); % for twist angle Shear modulus

%% Sea Level Analysis
rho1 = 1.225;                      % kg/m^3     at sea level  
nz = 100;                         % number of  ncrements along z-axis

% PHAA 
n1 = 4.4;                          %            load factor of PHAA
V1 = 59.5350;                      % m/s        velocity at PHAA
AoA1 = (19/180)*pi;                % rad        from XFoil results
Cd = 0.00586;                     % drag coefficient at zero lift
CM = 0.0006;                      % pitch moment coefficient from XFoil

% PLAA 
n2 = 4.4;
V2 = 95.8;
AoA2 = 0.070364436;


% NHAA
n3 = -1.76;
V3 = 39.8;
AoA3 = (-19/180)*pi;

% Dive Gust
n4 = -1.1106;
V4 = 95.8;
AoA4 = -0.065867042;

% Cruise Gust 
n5 =  -1.8142;
V5 = 63.9;
AoA5 =-0.139218653;



%% Altitude Analysis

rho6 = 0.782043;              % kg/m^3     at 14,600  

% PHAA 
n6 = 4.4;                     %  load factor of PHAA
V6 = 76;                      % m/s        velocity at PHAA
AoA6 = (18/180)*pi;           % rad        from XFoil results

% PLAA 
n7 = 4.4;
V7 = 95.8;
AoA7 = 0.131935904;

% NHAA
n8 = -1.76;
V8 = 51.3;
AoA8 = (-18/180)*pi;

% Dive Gust
n9 = -1.3279;
V9 = 95.8;
AoA9 = -0.089625889;
 
% Cruise Gust 
n10 =  -2.1039;
V10 = 63.9;
AoA10 = -0.221178085;

%%
n = [n1 n2 n3 n4 n5 n6 n7 n8 n9 n10];
V = [V1 V2 V3 V4 V5 V6 V7 V8 V9 V10];
aoa = [AoA1 AoA2 AoA3 AoA4 AoA5 AoA6 AoA7 AoA8 AoA9 AoA10];

%% Combined Drag Distribution Graphs
figure(1)
hold on % SHOULD BE 5
for i=1:1
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho1,V(i),aoa(i),Cd,CM,nz);
plot(z,d,'linewidth',2)
end
xlabel('z (m)')
ylabel('Drag Distribution Sea Level (N/m)')
legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

% figure(2)
% hold on
% for i=6:10
% [wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho6,V(i),aoa(i),Cd,CM,nz);
% plot(z,d,'linewidth',2)
% end
% xlabel('z (m)')
% ylabel('Drag Distribution Altitude (N/m)')
% legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

%% Combined Lift Distribution Graphs 
figure(3)
hold on % SHOULD BE 5
for i=1:1
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho1,V(i),aoa(i),Cd,CM,nz);
plot(z,l,'linewidth',2)
end
xlabel('z (m)')
ylabel('Lift Distribution Sea Level (N/m)')
legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

% figure(4)
% hold on
% for i=6:10
% [wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho6,V(i),aoa(i),Cd,CM,nz);
% plot(z,l,'linewidth',2)
% end
% xlabel('z (m)')
% ylabel('Lift Distribution Altitude (N/m)')
% legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

%% Combined Wx and Wy 

figure(5)
hold on % SHOULD BE 5
for i=1:1
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho1,V(i),aoa(i),Cd,CM,nz);
plot(z,wx0,'linewidth',2)
end
xlabel('z (m)')
ylabel('X Load Distribution Sea Level (N/m)')
legend('Wx PHAA','Wx PLAA','Wx NHAA','Wx Negative Dive Gust','Wx Negative Cruise Gust');

figure(6)
hold on % SHOULD BE 5
for i=1:1
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho1,V(i),aoa(i),Cd,CM,nz);
plot(z,wy0,'linewidth',2)
end
xlabel('z (m)')
ylabel('Y Load Distribution Sea Level (N/m)')
legend('Wy PHAA','Wy PLAA','Wy NHAA','Wy Negative Dive Gust','Wy Negative Cruise Gust');


% figure(7)
% hold on
% for i=6:10
% [wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho6,V(i),aoa(i),Cd,CM,nz);
% plot(z,wx0,'linewidth',2)
% end
% xlabel('z (m)')
% ylabel('X Load Distribution Altitude (N/m)')
% legend('Wx PHAA','Wx PLAA','Wx NHAA','Wx Negative Dive Gust','Wx Negative Cruise Gust');
% 
% figure(8)
% hold on
% for i=6:10
% [wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho6,V(i),aoa(i),Cd,CM,nz);
% plot(z,wy0,'linewidth',2)
% end
% xlabel('z (m)')
% ylabel('Y Load Distribution Altitude (N/m)')
% legend('Wy PHAA','Wy PLAA','Wy NHAA','Wy Negative Dive Gust','Wy Negative Cruise Gust');

%% airfoil section properties
sh = 0.08; %m
st = 0.025; %m spar thickness
kt = 0.001016; %m skin thickness

A_cap = 5e-5;
A_str = 2.5e-5;
t_spar = st;
t_skin = kt;
t_plate = t_skin;
% locations of spars, spar caps and stringers (nose at the origin of the coordinate)
x_spar0 = .5;                       % front spar (2 cell beam)
x_strU0 = [0 .15 .3 .7 .95];                     % upper surface
x_strL0 = [0 .15 .3 .7 .95];                     % lower surface
% new coordinate with origin at the centroid is used for the output below
[c,Ixx,Iyy,Ixy,x,yU,yL,x_strU,x_strL,x_boomU,x_boomL,L_boomU,L_boomL,x_spar,h_spar,i_spar,dx,y_strU,y_strL] = airfoil_section(A_cap,A_str,t_spar,t_skin,x_spar0,x_strU0,x_strL0);

%x_boomU = x_skinU = position of elements (nodes). Nodes at stringers and
%spar caps and in between each
%x_boomL = x_skinL =  same as above but for bottom side
%L_boomU = L_skinU = length of skin between nodes top
%L_boomL = L_skinL = length of skin between nodes bottom

%% calculate shear force, moments, and displacement along z-axis from root to tip
    if i>5
        rho = rho6;
    end

% Sea Level
for i = 1:1
  
  rho = rho1;
 
  [wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho,V(i),aoa(i),Cd,CM,nz);
      
    % Note wx0 and wy0 here are defined from root to tip.

wy = zeros(1,nz+1);
wx = zeros(1,nz+1);

for i = 1:nz+1
    wy(i) = wy0(nz+2-i);
    wx(i) = wx0(nz+2-i);
end
    
L = b/2;                                  % m   half span
z = 0:L/nz:L;

% Force y moment x
Sy = zeros(1,nz+1);
Sy_sum = zeros(1,nz+1);

Mx = zeros(1,nz+1);
Mx_sum = zeros(1,nz+1);

for i = 1:nz;
    Sy(i+1) = -(z(i+1)-z(i))*(wy(i+1) + wy(i))/2; % Force at each incriment
    Sy_sum(i+1) = Sy_sum(i) + Sy(i+1); % Total force up to incriment 
    Mx(i+1) = (z(i+1)-z(i))*(Sy_sum(i+1) + Sy_sum(i))/2; % Moment at each incriment
    Mx_sum(i+1) = Mx_sum(i) + Mx(i+1); % Total moment up to incriment
end

% z-->L-z: from root to tip
Mx0 = zeros(1,nz+1);
Sy0 = zeros(1,nz+1);
for i = 1:nz+1
    Mx0(i) = Mx_sum(nz+2-i);
    Sy0(i) = Sy_sum(nz+2-i);
end

%% Force x moment y
Sx = zeros(1,nz+1);
Sx_sum = zeros(1,nz+1);

My = zeros(1,nz+1);
My_sum = zeros(1,nz+1);

for i = 1:nz;
    Sx(i+1) = -(z(i+1)-z(i))*(wx(i+1) + wx(i))/2; % Force at each incriment
    Sx_sum(i+1) = Sx_sum(i) + Sx(i+1); % Total force up to incriment 
    My(i+1) = (z(i+1)-z(i))*(Sx_sum(i+1) + Sx_sum(i))/2; % Moment at each incriment
    My_sum(i+1) = My_sum(i) + My(i+1); % Total moment up to incriment
end

% z-->L-z: from root to tip
My0 = zeros(1,nz+1);
Sx0 = zeros(1,nz+1);
for i = 1:nz+1
    My0(i) = My_sum(nz+2-i);
    Sx0(i) = Sx_sum(nz+2-i);
end
Sx0 = -Sx0;
Sy0 = -Sy0;
%% Moment, Shear, and Distributed Load Graphs

figure(9)
hold on;
plot(z,Mx0./1e3,'Linewidth',2)
xlabel('z (m)')
ylabel('X Moments at Sea Level (kN*m)')
legend('Mx PHAA','Mx PLAA','Mx NHAA','Mx Negative Dive Gust','Mx Negative Cruise Gust');

figure(10)
hold on;
plot(z, My0./1e3, 'Linewidth',2)
xlabel('z (m)')
ylabel('Y Momemts at Sea Level (kN*m)')
legend('My PHAA','My PLAA','My NHAA','My Negative Dive Gust','My Negative Cruise Gust');


figure(11)
hold on;
plot(z,Sx0./1e3,'Linewidth',2)
xlabel('z (m)')
ylabel('X Shear Force at Sea Level (kN)')
legend('Sx PHAA','Sx PLAA','Sx NHAA','Sx Negative Dive Gust','Sx Negative Cruise Gust');


figure(12)
hold on;
plot(z,Sy0./1e3,'Linewidth',2)
xlabel('z (m)')
ylabel('Y Shear Force at Sea Level (kN)')
legend('Sy PHAA','Sy PLAA','Sy NHAA','Sy Negative Dive Gust','Sy Negative Cruise Gust');



%% calculate deflection u & v along z axis from root to tip
                           % modulus of the material
K = 1/(E*(Ixx*Iyy - Ixy^2));
ddu = -K.*(-Ixy.*Mx0 + Ixx.*My0);
ddv = -K.*(Iyy.*Mx0 - Ixy.*My0);

du = cumtrapz(z,ddu);
dv = cumtrapz(z,ddv);

u = cumtrapz(z,du);
v = cumtrapz(z,dv);

figure(13)
hold on;
plot(z,u,'Linewidth',2)
legend('u PHAA','u PLAA','u NHAA','u Negative Dive Gust','u Negative Cruise Gust');
xlabel('z (m)')
ylabel('X Deflection at Sea Level (m)')

figure(14)
hold on;
plot(z,v,'Linewidth',2)
legend('v PHAA','v PLAA','v NHAA','v Negative Dive Gust','v Negative Cruise Gust');
xlabel('z (m)')
ylabel('Y Deflection at Sea Level (m)')

% Stress_zz along contour of airfoil and spars 
nspars = 10; % points along spars 
yspar1 = yL(i_spar(1)):(h_spar(1))/nspars:yU(i_spar(1)); % points along spar 1
xspar1 = x(i_spar(1))*ones(1,nspars+1); % corresponding x position for spar 1
yspar2 = yL(i_spar(2)):(h_spar(2))/nspars:yU(i_spar(2)); % points along spar 2
xspar2 = x(i_spar(2))*ones(1,nspars+1); % corresponding x position for spar 1

xs = [flip(x),x,xspar1,xspar2]; % array of x coordinates going from te to te then 10 poitns for front and rear spar each
ys = [flip(yU),yL,yspar1,yspar2];
Sig_z = zeros(1,length(xs));

for i = 1:length(xs)
    
Sig_z(i) = (Mx0(1)*(Iyy*ys(i) - Ixy*xs(i)) + My0(1)*(Ixx*xs(i) - Ixy*ys(i)))/(Ixx*Iyy-Ixy^2);

end

figure(15)
hold on;
scatter(xs,Sig_z)
ylabel('Stress_{zz} at Sea Level(N/m^2)');
xlabel('x (m)');
legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

end

    
for i = 6:10
%          rho = rho6;
%  
%   [wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho,V(i),aoa(i),Cd,CM,nz);
%       
%     % Note wx0 and wy0 here are defined from root to tip.
% 
% wy = zeros(1,nz+1);
% wx = zeros(1,nz+1);
% for i = 1:nz+1
%     wy(i) = wy0(nz+2-i);
%     wx(i) = wx0(nz+2-i);
% end
%     
% L = b/2;                                  % m   half span
% z = 0:L/nz:L;
% 
% % Force y moment x
% Sy = zeros(1,nz+1);
% Sy_sum = zeros(1,nz+1);
% 
% Mx = zeros(1,nz+1);
% Mx_sum = zeros(1,nz+1);
% 
% for i = 1:nz;
%     Sy(i+1) = -(z(i+1)-z(i))*(wy(i+1) + wy(i))/2; % Force at each incriment
%     Sy_sum(i+1) = Sy_sum(i) + Sy(i+1); % Total force up to incriment 
%     Mx(i+1) = (z(i+1)-z(i))*(Sy_sum(i+1) + Sy_sum(i))/2; % Moment at each incriment
%     Mx_sum(i+1) = Mx_sum(i) + Mx(i+1); % Total moment up to incriment
% end
% 
% % z-->L-z: from root to tip
% Mx0 = zeros(1,nz+1);
% Sy0 = zeros(1,nz+1);
% for i = 1:nz+1
%     Mx0(i) = Mx_sum(nz+2-i);
%     Sy0(i) = Sy_sum(nz+2-i);
% end
% 
% %% Force x moment y
% Sx = zeros(1,nz+1);
% Sx_sum = zeros(1,nz+1);
% 
% My = zeros(1,nz+1);
% My_sum = zeros(1,nz+1);
% 
% for i = 1:nz;
%     Sx(i+1) = -(z(i+1)-z(i))*(wx(i+1) + wx(i))/2; % Force at each incriment
%     Sx_sum(i+1) = Sx_sum(i) + Sx(i+1); % Total force up to incriment 
%     My(i+1) = (z(i+1)-z(i))*(Sx_sum(i+1) + Sx_sum(i))/2; % Moment at each incriment
%     My_sum(i+1) = My_sum(i) + My(i+1); % Total moment up to incriment
% end
% 
% % z-->L-z: from root to tip
% My0 = zeros(1,nz+1);
% Sx0 = zeros(1,nz+1);
% for i = 1:nz+1
%     My0(i) = My_sum(nz+2-i);
%     Sx0(i) = Sx_sum(nz+2-i);
% end
% 
% %% Moment, Shear, and Distributed Load Graphs
% 
% figure(16)
% hold on;
% plot(z,Mx0./1e3,'Linewidth',2)
% legend('Mx PHAA','Mx PLAA','Mx NHAA','Mx Negative Dive Gust','Mx Negative Cruise Gust');
% xlabel('z (m)')
% ylabel(strcat('X Moments at Altitude (kN/m)'))
% 
% figure(17)
% hold on;
% plot(z,My0./1e3,'Linewidth',2)
% legend('My PHAA','My PLAA','My NHAA','My Negative Dive Gust','My Negative Cruise Gust');
% xlabel('z (m)')
% ylabel(strcat('Y Moments at Altitude (kN/m)'))
% 
% figure(18)
% hold on;
% plot(z,-Sx0./1e3,'Linewidth',2)
% xlabel('z (m)')
% ylabel(strcat('X Shear Force at Altitude (kN)'))
% legend('Sx PHAA','Sx PLAA','Sx NHAA','Sx Negative Dive Gust','Sx Negative Cruise Gust');
% 
% figure(19)
% hold on;
% plot(z,-Sy0./1e3,'Linewidth',2)
% xlabel('z (m)')
% ylabel(strcat('Y Shear Force at Altitude (kN)'))
% legend('Sy PHAA','Sy PLAA','Sy NHAA','Sy Negative Dive Gust','Sy Negative Cruise Gust');
% 
% 
% %% calculate deflection u & v along z axis from root to tip
% E = 70e9;                             % modulus of the material
% K = 1/(E*(Ixx*Iyy - Ixy^2));
% ddu = -K.*(-Ixy.*Mx0 + Ixx.*My0);
% ddv = -K.*(Iyy.*Mx0 - Ixy.*My0);
% 
% du = cumtrapz(z,ddu);
% dv = cumtrapz(z,ddv);
% 
% u = cumtrapz(z,du);
% v = cumtrapz(z,dv);
% 
% figure(20)
% hold on;
% plot(z,u,'Linewidth',2)
% legend('u PHAA','u PLAA','u NHAA','u Negative Dive Gust','u Negative Cruise Gust');
% xlabel('z (m)')
% ylabel('X Deflection at Altitude (m)')
% 
% figure(21)
% hold on;
% plot(z,v,'Linewidth',2)
% legend('v PHAA','v PLAA','v NHAA','v Negative Dive Gust','v Negative Cruise Gust');
% xlabel('z (m)')
% ylabel('Y Deflection at Altitude (m)')
% 
% 
% 
% Sig_z = (Mx0*(Iyy*wing.c_y - Ixy*wing.c_x) + My0*(Ixx*wing.c_x - Ixy*wing.c_y))/(Ixx*Iyy-Ixy^2);
% 
% figure (22)
% hold on;
% plot(z,Sig_z,'Linewidth',2)
% xlabel('z (m)')
% ylabel('Stress_{zz} at Altitude(N/m^2)');
% legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');
end
   
%% shear flow (BOOM AREAS)

nBU = length(x_boomU); % Number of nodes on upper side
i_BU = zeros(1,nBU);
i_BU(:) = round((x_boomU(:)-x_boomU(1))/dx)+1; % incriment location of upper boom in YU


% stress at each boom on the top surface at the root of the wing
sz_RBU = zeros(1,nBU);
sz_RBU(:) = Mx0(1)*(Iyy*yU(i_BU(:))-Ixy*x(i_BU(:)))/(Ixx*Iyy-Ixy^2) + My0(1)*(Ixx*x(i_BU(:))-Ixy*yU(i_BU(:)))/(Ixx*Iyy-Ixy^2);
if_strU = ismember(x_boomU,x_strU); % Array of 1 and 0 saying which in elements in x_boomU ( node locations) show up in x_strU (String locations)
if_capU = ismember(x_boomU,x_spar);  % Array of 1 and 0 saying which in elements in x_boomU ( node locations) show up in x_spar (String locations)

% Boom area upper portion
BU = zeros(1,nBU);
for i = 1:nBU % Places stringers and caps in boom array along proper section
    if if_strU(i) == 1
        BU(i) = A_str;
    elseif if_capU(i) == 1
        BU(i) = A_cap;
    end
end

% For first boom no area except for stringer or cap
for i = 2:nBU-1
    BU(i) = BU(i) + t_skin*L_boomU(i-1)/6*(2+sz_RBU(i-1)/sz_RBU(i)) + t_skin*L_boomU(i)/6*(2+sz_RBU(i+1)/sz_RBU(i)); % Stringer/Cap area + skin before and after
end
BU(end) = BU(end) + t_skin*L_boomU(end-1)/6*(2+sz_RBU(end-1)/sz_RBU(end)); % For last boom, stringer or cap area and only area left of boom 

nBL = length(x_boomL); % Number of nodes on lower side
i_BL = zeros(1,nBL);
i_BL(:) = round((x_boomL(:)-x_boomL(1))/dx)+1; % incriment location of lower boom in YL

% stress at each boom on the bottom surface at the root of the wing
sz_RBL = zeros(1,nBL);
sz_RBL(:) = Mx0(1)*(Iyy*yL(i_BL(:))-Ixy*x(i_BL(:)))/(Ixx*Iyy-Ixy^2) + My0(1)*(Ixx*x(i_BL(:))-Ixy*yL(i_BL(:)))/(Ixx*Iyy-Ixy^2);

if_strL = ismember(x_boomL,x_strL);  % Array of 1 and 0 saying which in elements in x_boomU ( node locations) show up in x_strU (String locations)
if_capL = ismember(x_boomL,x_spar); % Array of 1 and 0 saying which in elements in x_boomU ( node locations) show up in x_spar (spar locations)

% Boom area lower portion
BL = zeros(1,nBL); 
for i = 1:nBL  % Places stringers and caps in boom array along proper section (BL = array of node areas)
    if if_strL(i) == 1
        BL(i) = A_str;
    elseif if_capL(i) == 1
        BL(i) = A_cap;
    end
end

% First boom is just area of stringer 
for i = 2:nBL-1
    BL(i) = BL(i) + t_skin*L_boomL(i-1)/6*(2+sz_RBL(i-1)/sz_RBL(i)) + t_skin*L_boomL(i)/6*(2+sz_RBL(i+1)/sz_RBL(i));  % Stringer/Cap area + contribution due to skin before and after
end

% For last boom, stringer or cap area and only area left of boom 
BL(end) = BL(end) + t_skin*L_boomL(end-1)/6*(2+sz_RBL(end-1)/sz_RBL(end));

% the first stringer shared by upper and lower surface
BU(1) = BU(1) + t_skin*L_boomU(1)/6*(2+sz_RBU(2)/sz_RBU(1)) + t_skin*L_boomL(1)/6*(2+sz_RBL(2)/sz_RBL(1));
BL(1) = BU(1);

% add the contribution of spars
% index for spar location
i_BsparU = find(if_capU); % returns which elements of array are not zero (aka which node on top surface contains spars 
i_BsparL = find(if_capL); % Same as above but for lower surface 
for i = 1:2
    BU(i_BsparU(i)) = BU(i_BsparU(i)) + t_spar*h_spar(i)/6*(2+sz_RBL(i_BsparL(i))/sz_RBU(i_BsparU(i))); % at nodes with spars adds contribution due to spar (upper surface)
    BL(i_BsparL(i)) = BL(i_BsparL(i)) + t_spar*h_spar(i)/6*(2+sz_RBU(i_BsparU(i))/sz_RBL(i_BsparL(i))); % at nodes with spars adds contribution due to spar (lower surface) 
end

%% shear flow
% form cell 1: boom area, distance between booms and coordinates of each boom
B = [fliplr(BU),BL(2:end)];
x_boom = [fliplr(x_boomU),x_boomL(2:end)]; % array of x coordinate of booms from top TE to CCW to bottom TE
y_boom = [fliplr(yU(i_BU(:))),yL(i_BL(2:end))]; % array of y coordinate of booms from top TE to CCW to bottom TE
L_boom = [fliplr(L_boomU),L_boomL,h_spar(2)]; % array of skin lengths from top TE to lower TE and end spar

figure
plot(x_boom,y_boom,'ro','markersize',6)
hold on
plot(x,yU,'k',x,yL,'k','linewidth',1.5)
plot([x(end),x(end)],[yU(end),yL(end)],'b',[x_spar(1),x_spar(1)],[yU(i_spar(1)),yL(i_spar(1))],'b','linewidth',2)
ylim([-0.3 0.3]) 
xlabel('x (m)')
ylabel('y (m)')
legend('Boom!')
title('Boom Distribution')
grid on

% calculate area of triangle formed by two nodes on the airfoil profile and point(x_spar(1),0)
nq = length(L_boom);
A = zeros(1,nq);
for i = 1:nq-1
    A(i) = abs(x_boom(i)*(y_boom(i+1)-0) + x_boom(i+1)*(0-y_boom(i)) + x_spar(1)*(y_boom(i)-y_boom(i+1)))/2;
end
A(end) = abs(x_boom(end)*(y_boom(1)-0) + x_boom(1)*(0-y_boom(end)) + x_spar(1)*(y_boom(end)-y_boom(1)))/2;
Asum = sum(A); % Total area of cut section 

% calculate the area of each cell
i_A1 = find(ismember(x_boom,x_spar(1)));
A1 = A(i_A1(1):i_A1(2)-1); % Area sections from LE to spar
A1sum = sum(A1); % Area of cell closest to LE
A2sum = Asum - A1sum; % Area of cell nearest to TE

% calcualte qb at the root of the wing 
qb = zeros(1,nq);
% CHECK C1 and C2 
c1 = (Sy0(1)*Ixy - Sx0(1)*Ixx)/(Ixx*Iyy - Ixy^2);
c2 = (Sx0(1)*Ixy - Sy0(1)*Iyy)/(Ixx*Iyy - Ixy^2);
qb(1) = c1*B(1)*x_boom(1) + c2*B(1)*y_boom(1);
for i = 2:nq
   qb(i)= qb(i-1)+ c1*B(i)*x_boom(i) + c2*B(i)*y_boom(i);
end

% separate cell 1 and cell 2 
qb1 = qb(i_A1(1):i_A1(2)-1);
L_boom1 = L_boom(i_A1(1):i_A1(2)-1);

L1sum = sum(L_boom1);

qb2 = [qb(1:i_A1(1)-1),qb(i_A1(2):end)];
L_boom2 = [L_boom(1:i_A1(1)-1),L_boom(i_A1(2):end)];

% modify distance between booms at the rear spar to compensate the change of the thickness
L_boom2(end) = L_boom2(end)*t_skin/t_spar; 
L2sum = sum(L_boom2);

% equations
syms q01 q02

% Twisting
% Cell 1 sums
sumA1(1) = 0;
sumA1(2) = 0;
for i =  i_A1(1):i_A1(2)-1;
    sumA1(1) = sumA1(1) + B(i) * x_boom(i) * L_boom(i);
    sumA1(2) = sumA1(2) + B(i) * y_boom(i) * L_boom(i);
end

% Cell 2 sums 
sumA2(1) = 0;
sumA2(2) = 0;

for i = 1:i_A1(1)-1;
    sumA2(1) = sumA2(1) + B(i) * x_boom(i) * L_boom(i);
    sumA2(2) = sumA2(2) + B(i) * y_boom(i) * L_boom(i);
end 

for i = i_A1(2):length(L_boom);
    sumA2(1) = sumA2(1) + B(i) * x_boom(i) * L_boom(i);
    sumA2(2) = sumA2(2) + B(i) * y_boom(i) * L_boom(i);
end
    
% Moments
msum = 0;
for i = 1:nq-1;
    msum = msum + 2*qb(i)*A(i);
end

% Equating total moments
eq1 = -M0 + Sx0(1) * .25 * c + 2 * A1sum * q01 + 2 * A2sum * q02 + msum == 0;

% Equating angle of twist in both cells 

eq2 = (1/A1sum)*((q01/t_skin)*L1sum + (q01 - q02) * h_spar(1)/t_spar + c1*sumA1(1)/t_skin + c2*sumA1(2)/t_skin)...
    - (1/A2sum)*((q02/t_skin)*L2sum + (q02 - q01) * h_spar(2)/t_spar + c1*sumA2(1)/t_skin + c2*sumA2(2)/t_skin + q02*h_spar(2)/t_plate);

[q01,q02] = solve(eq1==0,eq2==0);

q01 = vpa(q01,4)
q02 = vpa(q02,4)


% shear flow along airfoil contour from top right corner to top right corner CCW 
q = zeros(1,nq);   
for i = 1:i_A1(1)-1 % top cell 2 
    q(i) = qb2(i) + q02;
end

for i = 1 : length(qb1) % cell 1 
    j = i + i_A1(1) - 1;
    q(j) = qb1(i) + q01;
end

for i = 1 : nq - i_A1(2) % bottom cell 2 
    j = i + i_A1(2)-1;
    k = j-length(qb2)+1;
    q(j) = qb2(k) + q02;
end

% Back plate/spar
q(end) = q02;


q_spar = q01 - q02;

dtdz = (.5/(G*A1sum))*((q01/t_skin)*L1sum + (q01 - q02) * h_spar(1)/t_spar + c1*sumA1(1)/t_skin + c2*sumA1(2)/t_skin);


% CHECKS

% CHECK 1: final element of qb goes to close to zero.... check
qbEnd = 100*qb(end)/max(abs(qb));

% CHECK 2,3: SX = sum of (qi*(x_i+1 - x_i) same for SY
SxShear = zeros(1,length(qb));
SyShear = zeros(1,length(qb)+1);

 for i = 1:length(SxShear)-1;
     SxShear(i) = q(i)*(x_boom(i+1) - x_boom(i));
     SyShear(i) = q(i)*(y_boom(i+1) - y_boom(i));
 end
 
 % ends that dont meet in x_boom and y_boom so do ends manuaully ..
 % SxCh assuming that end spars always mounted verticaly 
 SyShear(end-1) = q(end)*(y_boom(1)-y_boom(end));
 SyShear(end) = q_spar*(h_spar(1));
 
 SyCheck = [sum(SyShear) Sy0(1)];
 SxCheck = [sum(SxShear) Sx0(1)];
 
 % Percent difference between shear applied to estimated by shear flow 
 SyCheck(3) = 100*(SyCheck(1)-SyCheck(2))/(SyCheck(2));
 SxCheck(3) = 100*(SxCheck(1)-SxCheck(2))/(SxCheck(2));
 
 % CHECK 4: 
% shear stress tau
% please calculate the shear stress from the shear flow

shear = zeros(1,length(q));

for i = 1:length(shear);
    shear(i) = q(i)*kt*B(i);
end

TotalShearCheck = sum(shear);



[ribSpacing, maxStringerStress, strStress] = buckle(Mx0,My0,Ixx, Iyy, Ixy,x_strU,y_strU,x_strL,y_strL, A_str,E);