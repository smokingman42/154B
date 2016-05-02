clear all;
close all;
clc;

%structure dimensions
c = 1.5; %m
kt = 0.001016; %m skin thickness
sh = 0.08; %m
st = 0.0025; %m spar thickness
bl = 0.012; %m bracket height
bt = 0.0025; %m bracket thickness`
theta = 30; %degrees
alpha = 10; %degrees

wing = build_wing(c,bl,bt,sh,st,kt,theta);
wing = inertia_prop(wing,c,bl,bt,sh,st,kt,theta);
draw_wing(wing);
disp(wing.c_x)
disp(wing.c_y)
disp(wing.Ixx)
disp(wing.Iyy)
disp(wing.Ixy)

