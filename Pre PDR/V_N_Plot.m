function [] = V_N_Plot(W,c,S,V_C,V_D,C_Lmax,C_Lmin,C_La,rho,g,p)
 
%% Gust Loads
Ue_C = 50;                            %ft/s     %gust velocity at cruise
Ue_D = 25;                            %ft/s     %gust velocity at dive
mu = 2*(W/S)/(rho*c*C_La*g);
Kg = (0.88*mu)/(5.3+mu);
n_gust_C = 1 + ((Kg*C_La*Ue_C*V_C)/(498*(W/S))) ;
n_gust_D = 1 + ((Kg*C_La*Ue_D*V_D)/(498*(W/S))) ;
 
 
%% Convert to Metric
W = 4.44822162*W;                     % N
S = 0.09290304*S;                     % m^2
V_C = round(0.51444444*V_C,1);        % m/s
V_D = round(0.51444444*V_D,1);        % m/s
k = 0.1;                              % velocity increment 
V = 0:k:V_D;
rho = 515.378819*rho;                 % kg/m^3
 
%% Stall Speed
V_sp = sqrt((W*2)/(C_Lmax*rho*S));
V_sn = sqrt((W*2)/(abs(C_Lmin)*rho*S));
V_sp = round(V_sp,1);
V_sn = round(V_sn,1);
 
%% Positive Maneuvering Loads
n1 = 4.4;
PHAA = (((2*W*n1)/(rho*C_Lmax*S))^(1/2));   %Positive High Angle of Attack (Velocity at CL Max and n Max) 
PHAA = round(PHAA,1);

% Stall load factors
V1p = 0:k:PHAA;
n_pos1 = (0.5*rho*C_Lmax*V1p.^2*S)./W ;

%FAR Regulations
V2p = PHAA:k:V_D;
n_pos2 = 4.4*ones(1,length(V2p));

%For plot 
Vp=[V1p(1:length(V1p)-1),V2p];
n_pos = [n_pos1(1:length(n_pos1)-1),n_pos2];

% Negative Maneuvering Loads
n2 = -0.4*n1;
NHAA = sqrt((2*W*n2)/(rho*C_Lmin*S));

%Positive High Angle of Attack (Velcoity at CL Min and n min )
NHAA = round(NHAA,1);

%Stall load factor
V1n = 0:k:NHAA;
n_neg1 = (0.5*rho*C_Lmin*V1n.^2*S)./W;

%Far Regulations
V2n = NHAA:k:V_C;
n_neg2 = n2*ones(1,length(V2n));

%Remaining lower limit (n = -1 at V_D)
V3n = V_C:k:V_D;
n_neg3 = n2 + ((n2+1)/(V_C-V_D))*(V3n-V_C);

%Right bound
n_right = -1:.1:n1;
V_right = V_D*ones(length(n_right));

%For plot
Vn = [V1n(1:length(V1n)-1),V2n(1:length(V2n)-1),V3n];
n_neg = [n_neg1(1:length(n_neg1)-1),n_neg2(1:length(n_neg2)-1),n_neg3];

%% Positive Gust Loads
n_g1 = 1:(k*(n_gust_C-1)/(V_C)):n_gust_C; % (0,1) to (V_c, n_c)
n_g2 = n_gust_C:(k*(n_gust_D-n_gust_C)/(V_D-V_C)):n_gust_D; %(V_c,n_c) to (V_d, n_d)
n_gp1 = [n_g1(1:length(n_g1)-1),n_g2]; %Appends ^
n_gp2 = 1:(k*(n_gust_D-1)/(V_D)):n_gust_D; % (0,1) to (V_d, n_c)

%% velocity index
i_V_sp = round(V_sp/k) + 1; % Numer of elements in Velocity profile up to V stall Cl max +
i_V_sn = round(V_sn/k) + 1; % Numer of elements in Velocity profile up to V stall Cl max -
i_PHAA = round(PHAA/k) + 1; % Numer of elements in Velocity profile up to V stall Cl max +  n max
i_NHAA = round(NHAA/k) + 1; % Numer of elements in Velocity profile up to V stall Cl max - n min 
i_V_D = round(V_D/k) + 1; % Numer of elements in Velocity profile up to V dive 

%% Allowable Positive Envelope

% Stall loads
V1_allowp = V_sp:k:PHAA; % array of velocities from stall to PHAA 
n1_allowp = n_pos(i_V_sp:i_PHAA); % values of n from stall to PHAA

% Maneuver vs Gust Loads e(Positive)
V2_allowp = PHAA:k:V_D; % velocities from PHAA to V_D
n2_allowp = zeros(1,i_V_D-i_PHAA+1); % Corresponding n from PHAA to V_D (initialized)
 
% compares n from maneuvering and gust at those velocities and 
% greater of the two is saved for allowable 
for i = 1:length(n2_allowp)
    j = i - 1 + i_PHAA;
    if (n_pos(j) >= n_gp1(j))
        n2_allowp(i) = n_pos(j);
    else
        n2_allowp(i) = n_gp1(j);
    end
end

% For plot
V_allowp = [V1_allowp,V2_allowp];
n_allowp = [n1_allowp,n2_allowp];

% Allowable Negative Envelope

% Stall Loads 
V1_allown = V_sn:k:NHAA;
n1_allown = n_neg(i_V_sn:i_NHAA);

% Maneuver Loads 
V2_allown = NHAA:k:V_D;
n2_allown = zeros(1,i_V_D-i_NHAA+1);

% compares n from maneuvering and gust at those velocities and 
% greater of the two is saved for allowable 
for i = 1:length(n2_allown)
    j = i - 1 + i_NHAA;
    if(n_neg(j)<=-n_gp1(j)+2)
        n2_allown(i)=n_neg(j);
    else
        n2_allown(i)= -n_gp1(j)+2;
    end
end

% For plot
V_allown = [V1_allown,V2_allown];
n_allown = [n1_allown,n2_allown];


%% V-N Plot
figure('units','normalized','outerposition',[0 0 1 1]) ;
%Set the current figure to have 'Units' be 'normalized' and the 'outerposition', i.e. the 
%actual borders of the figure to be at the bottom left corner (0,0) and span the whole screen (1,1).
hold on
% Maneuvering Loads
p1 = plot(Vp,n_pos,'b','linewidth',5);
plot(Vn,n_neg,'b','linewidth',5);
plot(V_right,n_right,'b','linewidth',5);

% Gust Loads
p2 = plot(V,n_gp1,'r','linewidth',2);
plot(V,n_gp2,'r','linewidth',2);
plot(V,-n_gp1+2,'r','linewidth',2);
plot(V,-n_gp2+2,'r','linewidth',2);

% Acceptable Envelope
p3 = plot(V_allowp,n_allowp,'g','linewidth',4); % Positive load bound 
plot([V_sp,V_sp],[0,n_allowp(1)],'g','linewidth',4); % vertical line from 0 to 1 at V_sp 
plot([V_sp,V_sn],[0,0],'g','linewidth',4) % horizontal line from V_sp to V_Sn 
plot(V_allown,n_allown,'g','linewidth',4); % Negative load bound 
plot([V_sn,V_sn],[0,n_allown(1)],'g','linewidth',4); % Vertical line on V_sn from 0 to n(V_sn)
plot([V_D,V_D],[n_allown(length(n_allown)),n_allowp(length(n_allowp))],'g','linewidth',4); % Dive bound

% Critical Points
plot(PHAA, n1,'*','MarkerSize',15, 'MarkerEdgeColor', 'm','LineWidth',2);
text(PHAA+1, n1,'PHAA','VerticalAlignment','bottom','FontSize',14);
disp('PHAA')
disp(PHAA)
disp(n1)

plot(V_D, n1,'*','MarkerSize',15, 'MarkerEdgeColor', 'm','LineWidth',2);
text(V_D+1, n1,'PLAA','VerticalAlignment','bottom','FontSize',14);
disp('PLAA')
disp(V_D)
disp(n1)

plot(V_D, -1,'*','MarkerSize',15, 'MarkerEdgeColor', 'k','LineWidth',2);
text(V_D + 1, -1,'NLAA','VerticalAlignment','bottom','FontSize',14);

plot(NHAA, n2,'*','MarkerSize',15, 'MarkerEdgeColor', 'm','LineWidth',2);
text(NHAA+1, n2,'NHAA','VerticalAlignment','bottom','FontSize',14);
disp('NHAA')
disp(NHAA)
disp(n2)

plot(V_D, -n_gust_D+2,'*','MarkerSize',15, 'MarkerEdgeColor', 'm','LineWidth',2);
text(V_D+1, -n_gust_D+2,'Dive Gust','VerticalAlignment','top','FontSize',14);
disp('Dive Gust')
disp(V_D)
disp(-n_gust_D+2)


plot(V_C, -n_gust_C+2,'*','MarkerSize',15, 'MarkerEdgeColor', 'm','LineWidth',2);
text(V_C+1, -n_gust_C+2,'Cruise Gust','VerticalAlignment','top','FontSize',14);
disp('Cruise Gust')
disp(V_C)
disp(-n_gust_C+2)


plot(V_D, n_gust_D,'*','MarkerSize',15, 'MarkerEdgeColor', 'k','LineWidth',2);
text(V_D+1, n_gust_D,'Dive Gust','VerticalAlignment','top','FontSize',14);

plot(V_C, n_gust_C,'*','MarkerSize',15, 'MarkerEdgeColor', 'k','LineWidth',2);
text(V_C+1, n_gust_C+.1,'Cruise Gust','FontSize',14);


xlabel('Velocity, V [m/s]','fontsize',15)
ylabel('Load Factor, n [-]','fontsize',15)
set(gca,'fontsize',15)
h = legend([p1 p2 p3],'Maneuvering Loads','Gust Loads','Allowable Envelope','Location','northwest') ;
set(h,'fontsize',15);
grid on
switch p
    case 1
        title('V-N Diagram at Sea Level')
    case 2
        title('V-N Diagram at 14600 ft')
end
 
end
