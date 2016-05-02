%V-n Diagram

clear;                  %remove all variables from current workspace
close all;              %delete the current figure
clc;                    %clear command window
 
% Constants Metric
g_m = 9.81;               %m/s^2                 gravity
rho_sea_m = 1.225;        %kg/m^3             density of air at sea level
rho_alt_m = 0.782043;     %kg/m^3             density of air at 14,600 ft
e = 0.79 ;                %1                      Oswald efficiency

% Geometry Metric
b_m = 10.82;              %m                     span
c_m = 1.5;                %m                     chord length
S_m = b_m*c_m;            %m^2                   surface area
W_m = 1100;               %kg                    weight
AR= (b_m)^2/S_m;          %1                     aspect ratio

% Constants Imperial
g = convacc(g_m,'m/s^2','ft/s^2');                      %ft/s^2
rho_sea = convdensity(rho_sea_m,'kg/m^3','slug/ft^3');  %slug/ft^3 
rho_alt = convdensity(rho_alt_m,'kg/m^3','slug/ft^3');  %slug/ft^3

% Geometry Imperial
b = convlength(b_m,'m','ft');     %ft           span
c = convlength(c_m,'m','ft');     %ft           chord
S = b*c;                          %ft^2         surface area
W = convmass(W_m,'kg','lbm');     %lb           weight

% Given
V_C = 124.19;           %kts                    Cruise speed
V_D = 1.5*V_C;          %kts                    Dive speed
   
% 2D Lift Curve Slope: from XFoil
C_la_sea = 6.52937104;  %1/rad                  
C_la_alt = 6.5393104;   %1/rad                  
 
% 3D Lift Curve Slope
C_La_sea = C_la_sea/(1 + C_la_sea/(AR*pi*e));   %1/rad
C_La_alt = C_la_alt/(1 + C_la_alt/(AR*pi*e));   %1/rad

% 2D Lift alpha = 0
C_l0_sea = 0.2508;
C_l0_alt = 0.2501;


% 3d Lift alpha = 0
C_L0_sea =(C_La_sea/C_la_sea)*C_l0_sea;
C_L0_alt =(C_La_sea/C_la_sea)*C_l0_alt;
 
% 2D Lift Coefficient: from XFoil
C_lmax_sea = 1.8384;    %1                      
C_lmin_sea = -1.644;    %1
C_lmax_alt = 1.77;      %1
C_lmin_alt = -1.5551;    %1
 
% 3D Lift Coefficient
C_Lmax_sea = (C_La_sea/C_la_sea)*C_lmax_sea;    %1
C_Lmin_sea = (C_La_sea/C_la_sea )*C_lmin_sea;    %1
C_Lmax_alt = (C_La_alt/C_la_alt)*C_lmax_alt;    %1
C_Lmin_alt = (C_La_alt/C_la_alt)*C_lmin_alt;    %1
 
%% V-N Plots
disp('Sea Level')
% Sea Level
V_N_Plot(W,c,S,V_C,V_D,C_Lmax_sea,C_Lmin_sea,C_La_sea,rho_sea,g,1)
disp('Altitude')
% Altitude
V_N_Plot(W,c,S,V_C,V_D,C_Lmax_alt,C_Lmin_alt,C_La_alt,rho_alt,g,2)
