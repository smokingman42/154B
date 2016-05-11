function [c,Ixx,Iyy,Ixy,x,yU,yL,x_strU,x_strL,x_skinU,x_skinL,L_skinU,L_skinL,x_spar,h_spar,i_spar,dx] = airfoil_section(A_cap,A_str,t_spar,t_skin,x_spar,x_strU,x_strL)
%% airfoil section profile
% NACA 2412
m = 0.02;
p = 0.4;
t = 0.12;

c = 1.3462;             % chord length
nx = 500;               % number of increments
dx = c/nx;  
x = 0:dx:c;             % even spacing

yc = zeros(1,nx+1);
yt = zeros(1,nx+1);
yU0 = zeros(1,nx+1);
yL0 = zeros(1,nx+1);
theta = zeros(1,nx+1);
xb = p*c;
i_xb = xb/dx + 1;  

for i = 1:nx+1
    yt(i) = 5*t*c*(0.2969*sqrt(x(i)/c) - 0.1260*(x(i)/c) - 0.3516*(x(i)/c)^2 + 0.2843*(x(i)/c)^3 - 0.1015*(x(i)/c)^4);
    if i <= i_xb
        yc(i) = m*x(i)/p^2*(2*p - x(i)/c);
        theta(i) = atan(2*m/p^2*(p - x(i)/c));
    else
        yc(i) = m*(c - x(i))/(1-p)^2*(1 + x(i)/c - 2*p);
        theta(i) = atan(2*m/(1-p)^2*(p-x(i)/c));
    end
    yU0(i) = yc(i) + yt(i)*cos(theta(i));
    yL0(i) = yc(i) - yt(i)*cos(theta(i));
end

% The last 20% of the chord length of the airfoil was neglected under the assumption that
% this section contained flaps and ailerons, and would therefore not support aerodynamic loads.

i_xend = round(0.8*c/dx)+1;
x = x(1:i_xend);
yU = yU0(1:i_xend);
yL = yL0(1:i_xend);
% Here the airfoil profile is approximated by assuming xU0=x & xL0=x.

%% adding stringers, spar caps and spars
% spars & spar caps
x_spar = [x_spar,x(end)];         % add the rear spar
n_spar = length(x_spar);
i_spar = round(x_spar./dx)+1;        
h_spar = yU(i_spar) - yL(i_spar);
Cy_spar = (yU(i_spar) + yL(i_spar))/2;
A_spar = t_spar.*h_spar;

% stringers
n_strU = length(x_strU);
i_strU = round(x_strU./dx)+1;
n_strL = length(x_strL);
i_strL = round(x_strL./dx)+1;

% skins
% nodes include spar caps and stringers
x_nodeU = [x_spar,x_strU];
x_nodeU = sort(x_nodeU);

% skin adding one more node in between
x_skinU = zeros(1,2*length(x_nodeU)-1);
for i = 1:length(x_nodeU)-1
    x_skinU(2*i-1) = x_nodeU(i);
    x_skinU(2*i) = (x_nodeU(i) + x_nodeU(i+1))/2;
end
x_skinU(end) = x_nodeU(end);

i_skinU = round(x_skinU/dx)+1;
n_skinU = length(x_skinU)-1;
L_skinU = zeros(1,n_skinU);
A_skinU = zeros(1,n_skinU);
Cx_skinU = zeros(1,n_skinU);
Cy_skinU = zeros(1,n_skinU);
for i = 1:n_skinU
    L_skinU(i) = sqrt((yU(i_skinU(i+1)) - yU(i_skinU(i)))^2 + (x_skinU(i+1) - x_skinU(i))^2);
    A_skinU(i) = t_skin*L_skinU(i);
    Cx_skinU(i) = (x_skinU(i+1) + x_skinU(i))/2;
    Cy_skinU(i) = (yU(i_skinU(i+1)) + yU(i_skinU(i)))/2;
end

% lower part
x_nodeL = [x_spar,x_strL];
x_nodeL = sort(x_nodeL);

x_skinL = zeros(1,2*length(x_nodeL)-1);
for i = 1:length(x_nodeL)-1
    x_skinL(2*i-1) = x_nodeL(i);
    x_skinL(2*i) = (x_nodeL(i) + x_nodeL(i+1))/2;
end
x_skinL(end) = x_nodeL(end);

i_skinL = round(x_skinL/dx)+1;
n_skinL = length(x_skinL)-1;
L_skinL = zeros(1,n_skinL);
A_skinL = zeros(1,n_skinL);
Cx_skinL = zeros(1,n_skinL);
Cy_skinL = zeros(1,n_skinL);
for i = 1:n_skinL
    L_skinL(i) = sqrt((yL(i_skinL(i+1)) - yL(i_skinL(i)))^2 + (x_skinL(i+1) - x_skinL(i))^2);
    A_skinL(i) = t_skin*L_skinL(i);
    Cx_skinL(i) = (x_skinL(i+1) + x_skinL(i))/2;
    Cy_skinL(i) = (yL(i_skinL(i+1)) + yL(i_skinL(i)))/2;
end

%% centroid of the wing section
% initial value
Cx_sum = 0;
Cy_sum = 0;
A_sum = 0;

% Spars
for i = 1:n_spar
    Cx_sum = Cx_sum + x_spar(i)*A_spar(i);
    Cy_sum = Cy_sum + Cy_spar(i)*A_spar(i);
    A_sum = A_sum + A_spar(i);
end

% Upper Stringers
for i = 1:n_strU
    Cx_sum = Cx_sum + x_strU(i)*A_str;
    Cy_sum = Cy_sum + yU(i_strU(i))*A_str;
    A_sum = A_sum + A_str;
end

% Lower Stringers 
for i = 1:n_strL
    Cx_sum = Cx_sum + x_strU(i)*A_str;
    Cy_sum = Cy_sum + yL(i_strL(i))*A_str;
    A_sum = A_sum + A_str;
end

% Upper Spar Caps
for i = 1:n_spar
    Cx_sum = Cx_sum + x_spar(i)*A_cap;
    Cy_sum = Cy_sum + yU(i_spar(i))*A_cap;
    A_sum = A_sum + A_cap;
end

% Lower Spar Caps
for i = 1:n_spar
    Cx_sum = Cx_sum + x_spar(i)*A_cap;
    Cy_sum = Cy_sum + yL(i_spar(i))*A_cap;
    A_sum = A_sum + A_cap;
end

% Upper Skins
for i = 1:n_skinU
    Cx_sum = Cx_sum + Cx_skinU(i)*A_skinU(i);
    Cy_sum = Cy_sum + Cy_skinU(i)*A_skinU(i);
    A_sum = A_sum + A_skinU(i);
end

% Lower Skins
for i = 1:n_skinL
    Cx_sum = Cx_sum + Cx_skinL(i)*A_skinL(i);
    Cy_sum = Cy_sum + Cy_skinL(i)*A_skinL(i);
    A_sum = A_sum + A_skinL(i);
end


Cx = Cx_sum/A_sum;
Cy = Cy_sum/A_sum;

figure
plot(x,yU,'k',x,yL,'k','Linewidth',2);
ylim([-0.3 0.3])
hold on
plot(x_strU,yU(i_strU),'or',x_strL,yL(i_strL),'or','markersize',5);
plot([x_spar(1),x_spar(1)],[yU(i_spar(1)),yL(i_spar(1))],'b',[x(end),x(end)],[yU(end),yL(end)],'b','Linewidth',3);
plot(x_spar,yU(i_spar),'sg',x_spar,yL(i_spar),'sg','markersize',7);
plot(Cx,Cy,'m*')
ylabel('y (m)')
xlabel('x (m)')
title('Origin: Leading Edge');
grid on

%% Area moments of inertia
% Initialize
Ixx = 0;
Iyy = 0;
Ixy = 0;

% Spars
for i = 1:n_spar
    Ixx = Ixx + t_spar*h_spar(i)^3/12 + A_spar(i)*(Cy_spar(i)-Cy)^2;
    Iyy = Iyy + t_spar^3*h_spar(i)/12 + A_spar(i)*(x_spar(i)-Cx)^2;
    Ixy = Ixy + A_spar(i)*(Cy_spar(i)-Cy)*(x_spar(i)-Cx);
end

% Upper Spar Caps 
for i = 1:n_spar
    Ixx = Ixx + A_cap*(yU(i_spar(i)) - Cy)^2;
    Iyy = Iyy + A_cap*(x_spar(i) - Cx)^2;
    Ixy = Ixy + A_cap*(x_spar(i) - Cx)*(yU(i_spar(i)) - Cy);
end

% Lower Spar Caps 
for i = 1:n_spar
    Ixx = Ixx + A_cap*(yL(i_spar(i)) - Cy)^2;
    Iyy = Iyy + A_cap*(x_spar(i) - Cx)^2;
    Ixy = Ixy + A_cap*(x_spar(i) - Cx)*(yL(i_spar(i)) - Cy);
end

% Upper Stringers 
for i = 1:n_strU
    Ixx = Ixx + A_str*(yU(i_strU(i)) - Cy)^2;
    Iyy = Iyy + A_str*(x_strU(i) - Cx)^2;
    Ixy = Ixy + A_str*(yU(i_strU(i)) - Cy)*(x_strU(i) - Cx);
end

% Lower Stringers
for i = 1:n_strL
    Ixx = Ixx + A_str*(yL(i_strL(i)) - Cy)^2;
    Iyy = Iyy + A_str*(x_strL(i) - Cx)^2;
    Ixy = Ixy + A_str*(yL(i_strL(i)) - Cy)*(x_strL(i) - Cx);
end

% Upper Skin 
for i = 1:n_skinU
    Ixx = Ixx + A_skinU(i)*(Cy_skinU(i) - Cy)^2;
    Iyy = Iyy + A_skinU(i)*(Cx_skinU(i) - Cx)^2;
    Ixy = Ixy + A_skinU(i)*(Cy_skinU(i) - Cy)*(Cx_skinU(i) - Cx);
end

% Lower Skin 
for i = 1:n_skinL
    Ixx = Ixx + A_skinL(i)*(Cy_skinL(i) - Cy)^2;
    Iyy = Iyy + A_skinL(i)*(Cx_skinL(i) - Cx)^2;
    Ixy = Ixy + A_skinL(i)*(Cy_skinL(i) - Cy)*(Cx_skinL(i) - Cx);
end

%% Coordinate Transformation: New origin at centroid

% x Transformation 
for i =1:length(x)
    x(i) = x(i) - Cx;
end

% yU Transformation
for i = 1:length(yU)
    yU(i) = yU(i) - Cy;
end

% yL Transformation
for i = 1:length(yL)
    yL(i) = yL(i) - Cy;
end
    
% x_skinU Transformation
for i = 1:length(x_skinU)
    x_skinU(i) = x_skinU(i) - Cx;
end 

% x_skinL Transformation
for i = 1:length(x_skinL)
    x_skinL(i) = x_skinL(i) - Cx;
end 

% x_spar Transformation
for i = 1:length(x_spar)
    x_spar(i) = x_spar(i) - Cx;
end

% x_strU Transformation
for i = 1:length(x_strU)
    x_strU(i) = x_strU(i) - Cx;
end

% x_strL Transformation
for i = 1:length(x_strL)
    x_strL(i) = x_strL(i) - Cx;
end

Cx = 0;
Cy = 0;

figure
plot(x,yU,'k',x,yL,'k','Linewidth',2);
ylim([-0.3 0.3])
hold on
plot(x_strU,yU(i_strU),'or',x_strL,yL(i_strL),'or','markersize',5);
plot([x_spar(1),x_spar(1)],[yU(i_spar(1)),yL(i_spar(1))],'b',[x(end),x(end)],[yU(end),yL(end)],'b','Linewidth',3);
plot(x_spar,yU(i_spar),'sg',x_spar,yL(i_spar),'sg','markersize',7);
plot(Cx,Cy,'m*')
ylabel('y (m)')
xlabel('x (m)')
title('Origin: Centroid')
grid on

figure
plot(x,yU,'k',x,yL,'k','Linewidth',2);
ylim([-0.3 0.3])
hold on
plot(x_strU,yU(i_strU),'or',x_strL,yL(i_strL),'or','markersize',5);
plot([x_spar(1),x_spar(1)],[yU(i_spar(1)),yL(i_spar(1))],'b',[x(end),x(end)],[yU(end),yL(end)],'b','Linewidth',3);
plot(x_spar,yU(i_spar),'sg',x_spar,yL(i_spar),'sg','markersize',7);
plot(Cx,Cy,'m*')
ylabel('y (m)')
xlabel('x (m)')
grid on
scatter(x_skinU,yU(i_skinU)+.01,'v');
scatter(x_skinL,yL(i_skinL)-.01,'^');
title('Node Locations')


