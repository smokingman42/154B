%centroid and area of moment inertia
clear; clc; close;
%Bractet is considered as a point mass, which means the moments of inertia
%about axes through its own centroid are neglected.
%Skin and spars are assumed to have a rectanglar shape.
%3 spars, 8 brackets and 4 skin panels
%% input
x1 = 0;                                     %location of spar 1
x2 = 0.75;                                  %location of spar 2
x3 = 1.50;                                   %location of spar 3
t_spar = 0.0025;                            %spar thickness
t_skin = 0.001016;                          %skin thickness
A_bracket = 0.0025*0.015;                   %area of bracket
h = 0.08;                                   %height
h2 = 0.04;                                  %height of third spar
theta = 30/180*pi;                          %angle of trailing edge

L1 = x2-x1;                                 %skin length of leading edge
L2 = (x3-x2)/cos(theta);                    %skin length of trailing edge

%% Centroid
% centroid of each element
Cx_spar1 = x1;
Cx_spar2 = x2;
Cx_spar3 = x3;
Cy_spar1 = 0;
Cy_spar2 = 0;
Cy_spar3 = -(x3 - x2)*tan(theta);

Cx_skin1top = (x2 + x1)/2;
Cx_skin1bot = (x2 + x1)/2;
Cx_skin2top = (x3 + x2)/2;
Cx_skin2bot = (x3 + x2)/2;
Cy_skin1top = h/2;
Cy_skin1bot = -h/2;
Cy_skin2top = -(x3 - x2)/2*tan(theta) + h/2;
Cy_skin2bot = -(x3 - x2)/2*tan(theta) - h/2;

Cx_b1top = x1 + t_spar/2;
Cx_b1bot = x1 + t_spar/2;
Cy_b1top = h/2;
Cy_b1bot = -h/2;
Cx_b2top = x2 - t_spar/2;
Cx_b2bot = x2 - t_spar/2;
Cy_b2top = h/2;
Cy_b2bot = -h/2;
Cx_b3top = x2 + t_spar/2;
Cx_b3bot = x2 + t_spar/2;
Cy_b3top = h/2;
Cy_b3bot = -h/2;
Cx_b4top = x3 - t_spar/2;
Cx_b4bot = x3 - t_spar/2;
Cy_b4top = -(x3 - x2)*tan(theta) + h/2;
Cy_b4bot = -(x3 - x2)*tan(theta) - h/2;

%area of each element
A_spar = h*t_spar;
A_skin1 = L1*t_skin;
A_skin2 = L2*t_skin;

%centroid of the wing section
Cx = (Cx_spar1*A_spar + Cx_spar2*A_spar + Cx_spar3*A_spar...
     + A_skin1*Cx_skin1top + A_skin1*Cx_skin1bot...
     + A_skin2*Cx_skin2top + A_skin2*Cx_skin2bot...
     + A_bracket*Cx_b1top + A_bracket*Cx_b1bot...
     + A_bracket*Cx_b2top + A_bracket*Cx_b2bot...
     + A_bracket*Cx_b3top + A_bracket*Cx_b3bot...
     + A_bracket*Cx_b4top + A_bracket*Cx_b4bot)...
     /(3*A_spar + 2*A_skin1 + 2*A_skin2 + 8*A_bracket);

Cy = (Cy_spar1*A_spar + Cy_spar2*A_spar + Cy_spar3*A_spar...
     + A_skin1*Cy_skin1top + A_skin1*Cy_skin1bot...
     + A_skin2*Cy_skin2top + A_skin2*Cy_skin2bot...
     + A_bracket*Cy_b1top + A_bracket*Cy_b1bot...
     + A_bracket*Cy_b2top + A_bracket*Cy_b2bot...
     + A_bracket*Cy_b3top + A_bracket*Cy_b3bot...
     + A_bracket*Cy_b4top + A_bracket*Cy_b4bot)...
     /(3*A_spar + 2*A_skin1 + 2*A_skin2 + 8*A_bracket);
 
%% Area of moment inertia

Ixx_spar1 = t_spar*h^3/12 + A_spar*(Cy_spar1 - Cy)^2;
Iyy_spar1 = t_spar^3*h/12 + A_spar*(Cx_spar1 - Cx)^2;
Ixy_spar1 = A_spar*(Cx_spar1 - Cx)*(Cy_spar1 - Cy);
 
Ixx_spar2 = t_spar*h^3/12 + A_spar*(Cy_spar2 - Cy)^2;
Iyy_spar2 = t_spar^3*h/12 + A_spar*(Cx_spar2 - Cx)^2;
Ixy_spar2 = A_spar*(Cx_spar2 - Cx)*(Cy_spar2 - Cy);
 
Ixx_spar3 = t_spar*h^3/12 + A_spar*(Cy_spar3 - Cy)^2;
Iyy_spar3 = t_spar^3*h/12 + A_spar*(Cx_spar3 - Cx)^2;
Ixy_spar3 = A_spar*(Cx_spar3 - Cx)*(Cy_spar3 - Cy);
 
Ixx_b1top = A_bracket*(Cy_b1top - Cy)^2;
Iyy_b1top = A_bracket*(Cx_b1top - Cx)^2;
Ixy_b1top = A_bracket*(Cy_b1top - Cy)*(Cx_b1top - Cx);

Ixx_b1bot = A_bracket*(Cy_b1bot - Cy)^2;
Iyy_b1bot = A_bracket*(Cx_b1bot - Cx)^2;
Ixy_b1bot = A_bracket*(Cy_b1bot - Cy)*(Cx_b1bot - Cx);
 
Ixx_b2top = A_bracket*(Cy_b2top - Cy)^2;
Iyy_b2top = A_bracket*(Cx_b2top - Cx)^2;
Ixy_b2top = A_bracket*(Cy_b2top - Cy)*(Cx_b2top - Cx);

Ixx_b2bot = A_bracket*(Cy_b2bot - Cy)^2;
Iyy_b2bot = A_bracket*(Cx_b2bot - Cx)^2;
Ixy_b2bot = A_bracket*(Cy_b2bot - Cy)*(Cx_b2bot - Cx);

Ixx_b3top = A_bracket*(Cy_b3top - Cy)^2;
Iyy_b3top = A_bracket*(Cx_b3top - Cx)^2;
Ixy_b3top = A_bracket*(Cy_b3top - Cy)*(Cx_b3top - Cx);

Ixx_b3bot = A_bracket*(Cy_b3bot - Cy)^2;
Iyy_b3bot = A_bracket*(Cx_b3bot - Cx)^2;
Ixy_b3bot = A_bracket*(Cy_b3bot - Cy)*(Cx_b3bot - Cx);

Ixx_b4top = A_bracket*(Cy_b4top - Cy)^2;
Iyy_b4top = A_bracket*(Cx_b4top - Cx)^2;
Ixy_b4top = A_bracket*(Cy_b4top - Cy)*(Cx_b4top - Cx);

Ixx_b4bot = A_bracket*(Cy_b4bot - Cy)^2;
Iyy_b4bot = A_bracket*(Cx_b4bot - Cx)^2;
Ixy_b4bot = A_bracket*(Cy_b4bot - Cy)*(Cx_b4bot - Cx);

Ixx_skin1top = L1*t_skin^3/12 + A_skin1*(Cy_skin1top - Cy)^2;
Iyy_skin1top = L1^3*t_skin/12 + A_skin1*(Cx_skin1top - Cx)^2;
Ixy_skin1top = A_skin1*(Cx_skin1top - Cx)*(Cy_skin1top - Cy);
 
Ixx_skin1bot = L1*t_skin^3/12 + A_skin1*(Cy_skin1bot - Cy)^2;
Iyy_skin1bot = L1^3*t_skin/12 + A_skin1*(Cx_skin1bot - Cx)^2;
Ixy_skin1bot = A_skin1*(Cx_skin1bot - Cx)*(Cy_skin1bot - Cy);

Ixx_skin2top = L2^3*t_skin*(sin(theta))^2/12 + A_skin2*(Cy_skin2top - Cy)^2;
Iyy_skin2top = L2^3*t_skin*(cos(theta))^2/12 + A_skin2*(Cx_skin2top - Cx)^2;
Ixy_skin2top = -L2^3*t_skin*sin(2*theta)/24 + A_skin2*(Cx_skin2top - Cx)*(Cy_skin2top - Cy);
 
Ixx_skin2bot = L2^3*t_skin*(sin(theta))^2/12 + A_skin2*(Cy_skin2bot - Cy)^2;
Iyy_skin2bot = L2^3*t_skin*(cos(theta))^2/12 + A_skin2*(Cx_skin2bot - Cx)^2;
Ixy_skin2bot = -L2^3*t_skin*sin(2*theta)/24 + A_skin2*(Cx_skin2bot - Cx)*(Cy_skin2bot - Cy);

Ixx = Ixx_spar1 + Ixx_spar2 + Ixx_spar3...
    + Ixx_b1top + Ixx_b1bot + Ixx_b2top + Ixx_b2bot...
    + Ixx_b3top + Ixx_b3bot + Ixx_b4top + Ixx_b4bot...
    + Ixx_skin1top + Ixx_skin1bot + Ixx_skin2top + Ixx_skin2bot;

Iyy = Iyy_spar1 + Iyy_spar2 + Iyy_spar3...
    + Iyy_b1top + Iyy_b1bot + Iyy_b2top + Iyy_b2bot...
    + Iyy_b3top + Iyy_b3bot + Iyy_b4top + Iyy_b4bot...
    + Iyy_skin1top + Iyy_skin1bot + Iyy_skin2top + Iyy_skin2bot;

Ixy = Ixy_spar1 + Ixy_spar2 + Ixy_spar3...
    + Ixy_b1top + Ixy_b1bot + Ixy_b2top + Ixy_b2bot...
    + Ixy_b3top + Ixy_b3bot + Ixy_b4top + Ixy_b4bot...
    + Ixy_skin1top + Ixy_skin1bot + Ixy_skin2top + Ixy_skin2bot;





