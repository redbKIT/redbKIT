function [Disp] = rigid_motion(x, y, t, param)


switch param(1)
    case 1
        Disp = param(2) .*(x>0.1).*(x<0.9).*(y>0.1).*(y<0.9);
        
 
    case 2
        Disp = param(3) .*(x>0.1).*(x<0.9).*(y>0.1).*(y<0.9);
        
end


end