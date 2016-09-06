function [Disp] = deform_boundary(x, y, t, param)

cX = 0;
cY = 0;

theta =  cart2pol(x-cX,y-cY);

switch param(1)
    case 1
        Disp = (cX + cos(theta) + param(2)*cos(2*theta) - param(2) - x).*(x>-3).*(x<3).*(y>-3).*(y<3);
        
        %Disp = xx - x;

    case 2
        Disp = (cY + param(3) * sin(theta) - y).*(x>-3).*(x<3).*(y>-3).*(y<3);
        
        %Disp = yy - y;
end


end