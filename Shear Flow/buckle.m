function [ribSpacing,maxStress,Sig_zStr] = buckle(Mx0,My0, Ixx, Iyy, Ixy, x_strU, y_strU, x_strL, y_strL, A_str, E)

% Stringers modeled as squares
b = sqrt(A_str)
Ixx_str = (b^4)/12;

xstrLoc = [x_strU, x_strL];
ystrLoc = [y_strU, y_strL];
Sig_zStr = zeros(1,length(x_strU)+length(x_strL));

% Calc stresses at stringers 
for i = 1:length(xstrLoc) 
    Sig_zStr(i) = (Mx0(1)*(Iyy*ystrLoc(i) - Ixy*xstrLoc(i)) + My0(1)*(Ixx*xstrLoc(i) - Ixy*ystrLoc(i)))/(Ixx*Iyy-Ixy^2);
end

maxStress = abs(min(Sig_zStr));

ribSpacing = 2*pi*sqrt((E*Ixx_str)/(1.5*maxStress*A_str));


