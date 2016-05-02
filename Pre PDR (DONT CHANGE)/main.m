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

%% Sea Level Analysis
rho1 = 1.225;                      % kg/m^3     at sea level  
nz = 100;                         % number of increments along z-axis

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
hold on
for i=1:5
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho1,V(i),aoa(i),Cd,CM,nz);
plot(z,d,'linewidth',2)
end
xlabel('z (m)')
ylabel('Drag Distribution Sea Level (N/m)')
legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

figure(2)
hold on
for i=6:10
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho6,V(i),aoa(i),Cd,CM,nz);
plot(z,d,'linewidth',2)
end
xlabel('z (m)')
ylabel('Drag Distribution Altitude (N/m)')
legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

%% Combined Lift Distribution Graphs 
figure(3)
hold on
for i=1:5
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho1,V(i),aoa(i),Cd,CM,nz);
plot(z,l,'linewidth',2)
end
xlabel('z (m)')
ylabel('Lift Distribution Sea Level (N/m)')
legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

figure(4)
hold on
for i=6:10
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho6,V(i),aoa(i),Cd,CM,nz);
plot(z,l,'linewidth',2)
end
xlabel('z (m)')
ylabel('Lift Distribution Altitude (N/m)')
legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

%% Combined Wx and Wy 
figure(5)
hold on
for i=1:5
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho1,V(i),aoa(i),Cd,CM,nz);
plot(z,wx0,'linewidth',2)
plot(z,wy0,'linewidth',2)
end
xlabel('z (m)')
ylabel('Load Distribution Sea Level (N/m)')
legend('Wx PHAA','Wy PHAA','Wx PLAA','Wy PLAA','Wx NHAA','Wy NHAA','Wx Negative Dive Gust','Wy Negative Dive Gust','Wx Negative Cruise Gust','Wy Negative Cruise Gust');

figure(6)
hold on
for i=6:10
[wx0,wy0,M0,d,l,z] = get_wxwy(n(i),rho6,V(i),aoa(i),Cd,CM,nz);
plot(z,wx0,'linewidth',2)
plot(z,wy0,'linewidth',2)
end
xlabel('z (m)')
ylabel('Load Distribution Altitude (N/m)')
legend('Wx PHAA','Wy PHAA','Wx PLAA','Wy PLAA','Wx NHAA','Wy NHAA','Wx Negative Dive Gust','Wy Negative Dive Gust','Wx Negative Cruise Gust','Wy Negative Cruise Gust');



%% calculate shear force, moments, and displacement along z-axis from root to tip

wing = build_wing(c,bl,bt,sh,st,kt,theta);
wing = inertia_prop(wing,c,bl,bt,sh,st,kt,theta);

Ixx = wing.Ixx;
Iyy = wing.Iyy;
Ixy = wing.Ixy;

    if i>5
        rho = rho6;
    end


for i=1:5
  
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

%% Moment, Shear, and Distributed Load Graphs

figure(7)
hold on;
plot(z,Mx0./1e3,'Linewidth',2)
plot(z,My0./1e3,'Linewidth',2)
xlabel('z (m)')
ylabel('Moments at Sea Level (kN/m)')
legend('Mx PHAA','My PHAA','Mx PLAA','My PLAA','Mx NHAA','My NHAA','Mx Negative Dive Gust','My Negative Dive Gust','Mx Negative Cruise Gust','My Negative Cruise Gust');

figure(8)
hold on;
plot(z,-Sx0./1e3,'Linewidth',2)
plot(z,-Sy0./1e3,'Linewidth',2)
xlabel('z (m)')
ylabel('Shear Force at Sea Level (kN)')
legend('Sx PHAA','Sy PHAA','Sx PLAA','Sy PLAA','Sx NHAA','Sy NHAA','Sx Negative Dive Gust','Sy Negative Dive Gust','Sx Negative Cruise Gust','Sy Negative Cruise Gust');



%% calculate deflection u & v along z axis from root to tip
E = 70e9;                             % modulus of the material
K = 1/(E*(Ixx*Iyy - Ixy^2));
ddu = -K.*(-Ixy.*Mx0 + Ixx.*My0);
ddv = -K.*(Iyy.*Mx0 - Ixy.*My0);

du = cumtrapz(z,ddu);
dv = cumtrapz(z,ddv);

u = cumtrapz(z,du);
v = cumtrapz(z,dv);

figure(9)
hold on;
plot(z,u,'Linewidth',2)
plot(z,v,'Linewidth',2)
legend('u PHAA','v PHAA','u PLAA','v PLAA','u NHAA','v NHAA','u Negative Dive Gust','v Negative Dive Gust','u Negative Cruise Gust','v Negative Cruise Gust');
xlabel('z (m)')
ylabel('Deflection at Sea Level (m)')

Sig_z = (Mx0*(Iyy*wing.c_y - Ixy*wing.c_x) + My0*(Ixx*wing.c_x - Ixy*wing.c_y))/(Ixx*Iyy-Ixy^2);

figure(10)
hold on;
plot(z,Sig_z,'Linewidth',2)
xlabel('z (m)')
ylabel('Stress_{zz} at Sea Level(N/m^2)');
legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');

end

    
    for i = 6:10
         rho = rho6;
 
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

%% Moment, Shear, and Distributed Load Graphs

figure(11)
hold on;
plot(z,Mx0./1e3,'Linewidth',2)
plot(z,My0./1e3,'Linewidth',2)
legend('Mx PHAA','My PHAA','Mx PLAA','My PLAA','Mx NHAA','My NHAA','Mx Negative Dive Gust','My Negative Dive Gust','Mx Negative Cruise Gust','My Negative Cruise Gust');
xlabel('z (m)')
ylabel(strcat('Moments at Altitude (kN/m)'))

figure(12)
hold on;
plot(z,-Sx0./1e3,'Linewidth',2)
plot(z,-Sy0./1e3,'Linewidth',2)
xlabel('z (m)')
ylabel(strcat('Shear Force at Altitude (kN)'))
legend('Sx PHAA','Sy PHAA','Sx PLAA','Sy PLAA','Sx NHAA','Sy NHAA','Sx Negative Dive Gust','Sy Negative Dive Gust','Sx Negative Cruise Gust','Sy Negative Cruise Gust');



%% calculate deflection u & v along z axis from root to tip
E = 70e9;                             % modulus of the material
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
plot(z,v,'Linewidth',2)
legend('u PHAA','v PHAA','u PLAA','v PLAA','u NHAA','v NHAA','u Negative Dive Gust','v Negative Dive Gust','u Negative Cruise Gust','v Negative Cruise Gust');
xlabel('z (m)')
ylabel('Deflection at Altitude (m)')

Sig_z = (Mx0*(Iyy*wing.c_y - Ixy*wing.c_x) + My0*(Ixx*wing.c_x - Ixy*wing.c_y))/(Ixx*Iyy-Ixy^2);

figure(14)
hold on;
plot(z,Sig_z,'Linewidth',2)
xlabel('z (m)')
ylabel('Stress_{zz} at Altitude(N/m^2)');
legend('PHAA','PLAA','NHAA','Negative Dive Gust','Negative Cruise Gust');
    end
    
   






% % %====================================================================================================================
% % %% shear flow
% % % boom area upper part
% % nBU = length(x_boomU);
% % i_BU = zeros(1,nBU);
% % i_BU(:) = round((x_boomU(:)-x_boomU(1))/dx)+1;
% % 
% % % stress at each boom on the top surface at the root of the wing
% % sz_RBU = zeros(1,nBU);
% % sz_RBU(:) = Mx0(1)*(Iyy*yU(i_BU(:))-Ixy*x(i_BU(:)))/(Ixx*Iyy-Ixy^2) + My0(1)*(Ixx*x(i_BU(:))-Ixy*yU(i_BU(:)))/(Ixx*Iyy-Ixy^2);
% % 
% % if_strU = ismember(x_boomU,x_strU);
% % if_capU = ismember(x_boomU,x_spar);
% % 
% % BU = zeros(1,nBU);
% % for i = 1:nBU
% %     if if_strU(i) == 1
% %         BU(i) = A_str;
% %     elseif if_capU(i) == 1
% %         BU(i) = A_cap;
% %     end
% % end
% % 
% % for i = 2:nBU-1
% %     BU(i) = BU(i) + t_skin*L_boomU(i-1)/6*(2+sz_RBU(i-1)/sz_RBU(i)) + t_skin*L_boomU(i)/6*(2+sz_RBU(i+1)/sz_RBU(i));
% % end
% % BU(end) = BU(end) + t_skin*L_boomU(end-1)/6*(2+sz_RBU(end-1)/sz_RBU(end));
% % 
% % % boom area lower part
% % nBL = length(x_boomL);
% % i_BL = zeros(1,nBL);
% % i_BL(:) = round((x_boomL(:)-x_boomL(1))/dx)+1;
% % 
% % % stress at each boom on the bottom surface at the root of the wing
% % sz_RBL = zeros(1,nBL);
% % sz_RBL(:) = Mx0(1)*(Iyy*yL(i_BL(:))-Ixy*x(i_BL(:)))/(Ixx*Iyy-Ixy^2) + My0(1)*(Ixx*x(i_BL(:))-Ixy*yL(i_BL(:)))/(Ixx*Iyy-Ixy^2);
% % 
% % if_strL = ismember(x_boomL,x_strL);
% % if_capL = ismember(x_boomL,x_spar);
% % 
% % BL = zeros(1,nBL);
% % for i = 1:nBL
% %     if if_strL(i) == 1
% %         BL(i) = A_str;
% %     elseif if_capL(i) == 1
% %         BL(i) = A_cap;
% %     end
% % end
% % 
% % for i = 2:nBL-1
% %     BL(i) = BL(i) + t_skin*L_boomL(i-1)/6*(2+sz_RBL(i-1)/sz_RBL(i)) + t_skin*L_boomL(i)/6*(2+sz_RBL(i+1)/sz_RBL(i));
% % end
% % 
% % BL(end) = BL(end) + t_skin*L_boomL(end-1)/6*(2+sz_RBL(end-1)/sz_RBL(end));
% % 
% % % the first stringer shared by upper and lower surface
% % BU(1) = BU(1) + t_skin*L_boomU(1)/6*(2+sz_RBU(2)/sz_RBU(1)) + t_skin*L_boomL(1)/6*(2+sz_RBL(2)/sz_RBL(1));
% % BL(1) = BU(1);
% % 
% % % add the contribution of spars
% % % index for spar location
% % i_BsparU = find(if_capU);
% % i_BsparL = find(if_capL);
% % for i = 1:2
% %     BU(i_BsparU(i)) = BU(i_BsparU(i)) + t_spar*h_spar(i)/6*(2+sz_RBL(i_BsparL(i))/sz_RBU(i_BsparU(i)));
% %     BL(i_BsparL(i)) = BL(i_BsparL(i)) + t_spar*h_spar(i)/6*(2+sz_RBU(i_BsparU(i))/sz_RBL(i_BsparL(i)));
% % end
% % 
% % %% shear flow
% % % form cell 1: boom area, distance between booms and coordinates of each boom
% % B = [fliplr(BU),BL(2:end)];
% % x_boom = [fliplr(x_boomU),x_boomL(2:end)];
% % y_boom = [fliplr(yU(i_BU(:))),yL(i_BL(2:end))];
% % L_boom = [fliplr(L_boomU),L_boomL,h_spar(2)];
% % 
% % figure
% % plot(x_boom,y_boom,'ro','markersize',6)
% % hold on
% % plot(x,yU,'k',x,yL,'k','linewidth',1.5)
% % plot([x(end),x(end)],[yU(end),yL(end)],'b',[x_spar(1),x_spar(1)],[yU(i_spar(1)),yL(i_spar(1))],'b','linewidth',2)
% % ylim([-0.3 0.3]) 
% % xlabel('x (m)')
% % ylabel('y (m)')
% % title('Boom Distribution')
% % grid on
% % 
% % % calculate area of triangle formed by two nodes on the airfoil profile and point(x_spar(1),0)
% % nq = length(L_boom);
% % A = zeros(1,nq);
% % for i = 1:nq-1
% %     A(i) = abs(x_boom(i)*(y_boom(i+1)-0) + x_boom(i+1)*(0-y_boom(i)) + x_spar(1)*(y_boom(i)-y_boom(i+1)))/2;
% % end
% % A(end) = abs(x_boom(end)*(y_boom(1)-0) + x_boom(1)*(0-y_boom(end)) + x_spar(1)*(y_boom(end)-y_boom(1)))/2;
% % Asum = sum(A);
% % 
% % % calculate the area of each cell
% % i_A1 = find(ismember(x_boom,x_spar(1)));
% % A1 = A(i_A1(1):i_A1(2)-1);
% % A1sum = sum(A1);
% % A2sum = Asum - A1sum;
% % 
% % % calcualte qb at the root of the wing 
% % %
% % %
% % % please finish this part
% % %
% % %
% % 
% % % separate cell 1 and cell 2 
% % qb1 = qb(i_A1(1):i_A1(2)-1);
% % L_boom1 = L_boom(i_A1(1):i_A1(2)-1);
% % L1sum = sum(L_boom1);
% % 
% % qb2 = [qb(1:i_A1(1)-1),qb(i_A1(2):end)];
% % L_boom2 = [L_boom(1:i_A1(1)-1),L_boom(i_A1(2):end)];
% % % modify distance between booms at the rear spar to compensate the change of the thickness
% % L_boom2(end) = L_boom2(end)*t_skin/t_spar; 
% % L2sum = sum(L_boom2);
% % 
% % % equations
% % syms q01 q02
% % % eq1: equating moments of applied shear and pitch moment to moments of internal shear flow
% % %
% % % please set up equation 1 here by using eq1=...
% % %
% % 
% % % eq2: set the angle of twist of each cell to be the same
% % %
% % % please set up equation 2 here by using eq2=...
% % %
% %   
% % [q01,q02] = solve(eq1==0,eq2==0);
% % 
% % 
% % % shear flow along airfoil contour from top right corner to top right corner CCW 
% % q = zeros(1,nq);   
% % for i = 1:i_A1(1)-1
% %     q(i) = qb2(i) + q02;
% % end
% % 
% % for i = 1 : length(qb1)
% %     j = i + i_A1(1) - 1;
% %     q(j) = qb1(i) + q01;
% % end
% % 
% % for i = 1 : nq - i_A1(2) + 1
% %     j = i + i_A1(2) - 1;
% %     k = i + i_A1(1) - 1;
% %     q(j) = qb2(k) + q02;
% % end
% % 
% % %
% % % please calculate the shear flow in the central spar here
% % %
% % 
% % % verification:
% % %
% % % please verify your results in 4 ways (see instructions)
% % %
% % 
% % 
% % % shear stress tau
% % %
% % % please calculate the shear stress from the shear flow
% % %
% %     
% %     
% %     
% %     