function [] = draw_wing(wing)
    
    hold on
    
    %draw brackets
    for i = 1:16
       bracket_x = [wing.bracket(i).x(1) wing.bracket(i).x(3) ...
           wing.bracket(i).x(4) wing.bracket(i).x(2) ...
           wing.bracket(i).x(1)];
       bracket_y = [wing.bracket(i).y(1) wing.bracket(i).y(3) ...
           wing.bracket(i).y(4) wing.bracket(i).y(2) ...
           wing.bracket(i).y(1)];
       
       plot(wing.bracket(i).c_x,wing.bracket(i).c_y,'x', 'MarkerEdgeColor', 'r','MarkerSize', 5)
       plot(bracket_x,bracket_y,'k');
    end
    
    %draw spars
    for i = 1:3
       spar_x = [wing.spar(i).x(1) wing.spar(i).x(3) wing.spar(i).x(4) ...
           wing.spar(i).x(2) wing.spar(i).x(1)];
       spar_y = [wing.spar(i).y(1) wing.spar(i).y(3) wing.spar(i).y(4) ...
           wing.spar(i).y(2) wing.spar(i).y(1)];
       
       plot(wing.spar(i).c_x,wing.spar(i).c_y,'x', 'MarkerEdgeColor', 'r','MarkerSize', 5)
       plot(spar_x,spar_y,'g');
    end
    
    %draw skins
    for i = 1:4
       skin_x = [wing.skin(i).x(1) wing.skin(i).x(3) wing.skin(i).x(4) ...
           wing.skin(i).x(2) wing.skin(i).x(1)];
       skin_y = [wing.skin(i).y(1) wing.skin(i).y(3) wing.skin(i).y(4) ...
           wing.skin(i).y(2) wing.skin(i).y(1)];
       
       plot(wing.skin(i).c_x,wing.skin(i).c_y,'x', 'MarkerEdgeColor', 'r','MarkerSize', 5)
       plot(skin_x,skin_y,'b');
    end
    
    %total centroid
    plot(wing.c_x, wing.c_y, 'x', 'MarkerEdgeColor', 'r','MarkerSize', 10)
    
    %format plot
    axis equal;
    title('Wing Cross Section');
    xlabel('Length (m)');
    ylabel('Length (m)');
    
    hold off
    
end