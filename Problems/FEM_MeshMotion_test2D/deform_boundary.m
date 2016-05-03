function [Disp] = deform_boundary(x, y, t, param)

cX = 0.5;
cY = 0.5;

theta =  cart2pol(x-cX,y-cY);

switch param(1)
    case 1
        Disp = (cX + param(2) * cos(theta) - x).*(x>0.1).*(x<0.9).*(y>0.1).*(y<0.9);
        
        %Disp = xx - x;

    case 2
        Disp = (cY + param(3) * sin(theta) - y).*(x>0.1).*(x<0.9).*(y>0.1).*(y<0.9);
        
        %Disp = yy - y;
end


end