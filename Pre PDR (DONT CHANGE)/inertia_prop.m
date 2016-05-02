function [wing] = inertia_prop(wing,c,bl,bt,sh,st,kt,theta)
    % Uses simplified geometry of brackets, stringers, spars and skin to
    % find approximations for inertia and centroid
    
    %slanted section length
    X = (-c*cosd(theta)+sqrt(c^2*cosd(theta)^2+3*c^2))/2;
        
    % Areas
    bracketArea = bl*bt;
    skinArea13 = kt*c/2;
    skinArea24 = kt*X;
    sparArea = sh*st;
    totalArea = 16*bracketArea+2*skinArea13+2*skinArea24+3*sparArea; 
    
  %%  Centroids
    c_xSparSum = 0;
    c_ySparSum = 0;
    
    c_xBracketSum = 0;
    c_yBracketSum = 0;
    
    % Skin
    for i = 1:4 
        wing.skin(i).c_x = mean(wing.skin(i).x);
        wing.skin(i).c_y = mean(wing.skin(i).y);
    end
    
    % Spars
    for i = 1:3
        wing.spar(i).c_x = mean(wing.spar(i).x);
        wing.spar(i).c_y = mean(wing.spar(i).y);
        
        c_xSparSum = c_xSparSum + wing.spar(i).c_x;
        c_ySparSum = c_ySparSum + wing.spar(i).c_y;
    end
     
    % Brackets 
    for i = 1:16
        wing.bracket(i).c_x = mean(wing.bracket(i).x);
        wing.bracket(i).c_y = mean(wing.bracket(i).y);
        c_xBracketSum = c_xBracketSum + wing.bracket(i).c_x;
        c_yBracketSum = c_yBracketSum + wing.bracket(i).c_y;
    end
   
   wing.c_x = ((wing.skin(1).c_x + wing.skin(3).c_x)*skinArea13 + (wing.skin(2).c_x + wing.skin(4).c_x)*skinArea24...
         + c_xSparSum*sparArea + c_xBracketSum*bracketArea)/totalArea;
   
   wing.c_y = ((wing.skin(1).c_y + wing.skin(3).c_y)*skinArea13 + (wing.skin(2).c_y + wing.skin(4).c_y)*skinArea24...
        + c_ySparSum*sparArea + c_yBracketSum*bracketArea)/totalArea;
   
  %% Ixx 
   
   % Skin
     % Non-Tilted
   wing.skin(1).Ixx = c/2*kt^3/12 + skinArea13*(wing.skin(1).c_y - wing.c_y)^2;     
   wing.skin(3).Ixx = c/2*kt^3/12 + skinArea13*(wing.skin(3).c_y - wing.c_y)^2;
    
     % Titlted
  wing.skin(2).Ixx = (X^3*kt*(sin(-theta))^2)/12 + skinArea24*(wing.skin(2).c_y - wing.c_y)^2;
  wing.skin(4).Ixx = (X^3*kt*(sin(-theta))^2)/12 + skinArea24*(wing.skin(4).c_y - wing.c_y)^2;
  
  % Spars
    % Non-Tilted
  wing.spar(1).Ixx = (st*sh^3/12) + sparArea*(wing.spar(1).c_y - wing.c_y)^2;
  wing.spar(2).Ixx = (st*sh^3/12) + sparArea*(wing.spar(2).c_y - wing.c_y)^2;
  
    %Tilted
  wing.spar(3).Ixx = st^3*sh*(sin(-theta))^2/12 + sparArea*(wing.spar(3).c_y - wing.c_y)^2;
  
  % Brackets 
    % Non-Tilted Brackets 1-4
    for i = 1:4
        wing.bracket(i).Ixx = (wing.bracket(i).x(3)-wing.bracket(i).x(1))*(wing.bracket(i).y(1)-wing.bracket(i).y(2))^3/12 + ...
            bracketArea*(wing.bracket(i).c_y - wing.c_y)^2;
    end
    
    % Non-Tilted Brackets 9-12
    for i = 9:12
         wing.bracket(i).Ixx = (wing.bracket(i).x(3)-wing.bracket(i).x(1))*(wing.bracket(i).y(1)-wing.bracket(i).y(2))^3/12 + ...
            bracketArea*(wing.bracket(i).c_y - wing.c_y)^2;
    end
    
    % Tilted Brackets  (vertical and horizontal tilted)
    for i = 1:4
        bracket_no_vertical = [5 7 13 15];
        bracket_no_horizontal = [6 8 14 16];
        
        j = bracket_no_vertical(i);
        k = bracket_no_horizontal(i);
        
        wing.bracket(j).Ixx = (bt^3)*bl*(sin(-theta))^2/12 + bracketArea*(wing.bracket(j).c_y - wing.c_y)^2;
        wing.bracket(k).Ixx = (bl^3)*bt*(sin(-theta))^2/12 + bracketArea*(wing.bracket(k).c_y - wing.c_y)^2;
    end
        
  %% Iyy
   
   % Skin
     % Non-Tilted
   wing.skin(1).Iyy = (c/2)^3*kt/12 + skinArea13*(wing.skin(1).c_x - wing.c_x)^2;     
   wing.skin(3).Iyy = (c/2)^3*kt/12 + skinArea13*(wing.skin(3).c_x - wing.c_x)^2;
    
     % Titlted
  wing.skin(2).Iyy = (X^3*kt*(cos(-theta))^2)/12 + skinArea24*(wing.skin(2).c_x - wing.c_x)^2;
  wing.skin(4).Iyy = (X^3*kt*(cos(-theta))^2)/12 + skinArea24*(wing.skin(4).c_x - wing.c_x)^2;
  
  % Spars
    % Non-Tilted
  wing.spar(1).Iyy = ((st^3)*sh/12) + sparArea*(wing.spar(1).c_x - wing.c_x)^2;
  wing.spar(2).Iyy = ((st^3)*sh/12) + sparArea*(wing.spar(2).c_x - wing.c_x)^2;
  
    %Tilted
  wing.spar(3).Iyy = st^3*sh*(cos(-theta))^2/12 + sparArea*(wing.spar(3).c_x - wing.c_x)^2;
  
  % Brackets 
    % Non-Tilted Brackets 1-4
    for i = 1:4
        wing.bracket(i).Iyy = (wing.bracket(i).x(3)-wing.bracket(i).x(1))^3*(wing.bracket(i).y(1)-wing.bracket(i).y(2))/12 + ...
            bracketArea*(wing.bracket(i).c_x - wing.c_x)^2;
    end
    
    % Non-Tilted Brackets 9-12
    for i = 9:12
         wing.bracket(i).Iyy = (wing.bracket(i).x(3)-wing.bracket(i).x(1))^3*(wing.bracket(i).y(1)-wing.bracket(i).y(2))/12 + ...
            bracketArea*(wing.bracket(i).c_y - wing.c_y)^2;
    end
    
    % Tilted Brackets  (vertical and horizontal tilted)
    for i = 1:4
        bracket_no_vertical = [5 7 13 15];
        bracket_no_horizontal = [6 8 14 16];
        
        j = bracket_no_vertical(i);
        k = bracket_no_horizontal(i);
        
        wing.bracket(j).Iyy = (bt^3)*bl*(cos(-theta))^2/12 + bracketArea*(wing.bracket(j).c_y - wing.c_y)^2;
        wing.bracket(k).Iyy = (bl^3)*bt*(cos(-theta))^2/12 + bracketArea*(wing.bracket(k).c_y - wing.c_y)^2;
    end
    
  %% Ixy
  
   % Skin 
     % Non-Tilted
  wing.skin(1).Ixy = skinArea13*(wing.skin(1).c_x - wing.c_x)*(wing.skin(1).c_y - wing.c_y);
  wing.skin(3).Ixy = skinArea13*(wing.skin(3).c_x - wing.c_x)*(wing.skin(3).c_y - wing.c_y);
        
     % Tilted
  wing.skin(2).Ixy = X^3*kt*(sin(-2*theta))/24 + skinArea24*(wing.skin(2).c_x - wing.c_x)*(wing.skin(2).c_y - wing.c_y);
  wing.skin(4).Ixy = X^3*kt*(sin(-2*theta))/24 + skinArea24*(wing.skin(4).c_x - wing.c_x)*(wing.skin(4).c_y - wing.c_y); 
    
  % Spars 
     % Non - Tilted
  wing.spar(1).Ixy = sparArea*(wing.spar(1).c_x - wing.c_x)*(wing.spar(1).c_y - wing.c_y);
  wing.spar(2).Ixy = sparArea*(wing.spar(2).c_x - wing.c_x)*(wing.spar(2).c_y - wing.c_y);
   
     % Tilted
  wing.spar(3).Ixy =(st^3)*sh*(sin(-2*theta))/24 + sparArea*(wing.spar(3).c_x - wing.c_x)*(wing.spar(3).c_y - wing.c_y);
    
  % Brackets
    % Non - Tilted 
    for i = 1:8
        bracket_no = [ 1 2 3 4 9 10 11 12];
        k = bracket_no(i);
        wing.bracket(k).Ixy = bracketArea*(wing.bracket(k).c_x - wing.c_x)*(wing.bracket(k).c_y - wing.c_y);
    end
    
    % Tilted
    for i = 1:4
        bracket_no_vertical = [5 7 13 15];
        bracket_no_horizontal = [6 8 14 16];
        
        j = bracket_no_vertical(i);
        k = bracket_no_horizontal(i);
        
        wing.bracket(j).Ixy = (bt^3)*bl*(sin(-2*theta))/24 + bracketArea*(wing.bracket(j).c_x - wing.c_x)*(wing.bracket(j).c_y - wing.c_y);
        wing.bracket(k).Ixy = (bl^3)*bt*(sin(-2*theta))/24 + bracketArea*(wing.bracket(k).c_x - wing.c_x)*(wing.bracket(k).c_y - wing.c_y);
    end
    
        
 %% Total inertia
 
  wing.Ixx = 0;
  wing.Iyy = 0;
  wing.Ixy = 0;
        
  % Skins
    for i = 1:4
        wing.Ixx = wing.Ixx + wing.skin(i).Ixx;
        wing.Iyy = wing.Iyy + wing.skin(i).Iyy;
        wing.Ixy = wing.Ixy + wing.skin(i).Ixy;
    end
    
  % Spars 
    for i = 1:3
        wing.Ixx = wing.Ixx + wing.spar(i).Ixx;
        wing.Iyy = wing.Ixx + wing.spar(i).Iyy;
        wing.Ixy = wing.Ixy + wing.spar(i).Ixy;
    end
    
  % Brackets 
    for i = 1:16
        wing.Ixx = wing.Ixx + wing.bracket(i).Ixx;
        wing.Iyy = wing.Iyy + wing.bracket(i).Iyy;
        wing.Ixy = wing.Ixy + wing.bracket(i).Ixy;
    end
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
     
    