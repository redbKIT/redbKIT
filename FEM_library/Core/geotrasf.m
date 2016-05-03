function [detjac, invjac, h] = geotrasf(dim, vertices, elements)
%GEOTRASF linear 2\3 dimensional geometrical transformation
%
%   F. Saleri 24-08-01, F. Negri 18.11.2014

noe    = size(elements,2);
detjac = zeros(1,noe);
invjac = zeros(noe, dim, dim);
h      = zeros(1,noe);

switch dim
    
    case 2
        
        % Corner point indices
        a1 = elements(1,:); a2 = elements(2,:); a3 = elements(3,:);
        
        % Triangle sides
        s13x = vertices(1,a1)-vertices(1,a3);
        s31y = vertices(2,a3)-vertices(2,a1);
        s32x = vertices(1,a3)-vertices(1,a2);
        s23y = vertices(2,a2)-vertices(2,a3);
        s21x = vertices(1,a2)-vertices(1,a1);
        s21y = vertices(2,a2)-vertices(2,a1);
        
        % Determinant of the Jacobian matrix with sign (two times the area with sign of T)
        detjac = s13x.*s23y-s31y.*s32x;
        
        uno_su_detjac = 1./detjac;
        
        % Jacobian elements
        dcdx = s31y.*uno_su_detjac;
        dcdy = s13x.*uno_su_detjac;
        dedx = -s21y.*uno_su_detjac;
        dedy = s21x.*uno_su_detjac;
        h    = max([s13x.^2+s31y.^2;s32x.^2+s23y.^2;s21x.^2+s21y.^2]);
        h    = sqrt(h);
        
        invjac        = zeros(size(elements,2), 2, 2);
        
        
        invjac(:,1,1) = dcdx;
        invjac(:,1,2) = dedx;
        invjac(:,2,1) = dcdy;
        invjac(:,2,2) = dedy;
        
        
        % Determinant of the Jacobian of the linear transformation
        detjac = abs(detjac);
        
    case 3
        
        % Vertices
        a1 = elements(1,:); a2 = elements(2,:); a3 = elements(3,:); a4 = elements(4,:);
        %
        % Lengths of the sides
        s21x = vertices(1,a2)-vertices(1,a1);
        s31x = vertices(1,a3)-vertices(1,a1);
        s41x = vertices(1,a4)-vertices(1,a1);
        s21y = vertices(2,a2)-vertices(2,a1);
        s31y = vertices(2,a3)-vertices(2,a1);
        s41y = vertices(2,a4)-vertices(2,a1);
        s21z = vertices(3,a2)-vertices(3,a1);
        s31z = vertices(3,a3)-vertices(3,a1);
        s41z = vertices(3,a4)-vertices(3,a1);
        pzkI   = s31y.*s41z-s41y.*s31z;
        pzkII  = s31x.*s41z-s41x.*s31z;
        pzkIII = s31x.*s41y-s41x.*s31y;
        %
        % Volumes (multiplied for 6)
        volumes = s21x.*pzkI-s21y.*pzkII+s21z.*pzkIII;
        
        uno_su_volumes = 1./volumes;
        %
        % Jacobian elements
        dcdx = pzkI.*uno_su_volumes;
        dcdy = -pzkII.*uno_su_volumes;
        dcdz = pzkIII.*uno_su_volumes;
        
        dedx = (s41y.*s21z-s21y.*s41z).*uno_su_volumes;
        dedy = (s21x.*s41z-s41x.*s21z).*uno_su_volumes;
        dedz = (s21y.*s41x-s21x.*s41y).*uno_su_volumes;
        dtdx = (s21y.*s31z-s31y.*s21z).*uno_su_volumes;
        dtdy = (s21z.*s31x-s21x.*s31z).*uno_su_volumes;
        dtdz = (s21x.*s31y-s31x.*s21y).*uno_su_volumes;
        
        s23x = vertices(1,a2)-vertices(1,a3);
        s43x = vertices(1,a4)-vertices(1,a3);
        s24x = vertices(1,a2)-vertices(1,a4);
        s23y = vertices(2,a2)-vertices(2,a3);
        s43y = vertices(2,a4)-vertices(2,a3);
        s24y = vertices(2,a2)-vertices(2,a4);
        s23z = vertices(3,a2)-vertices(3,a3);
        s43z = vertices(3,a4)-vertices(3,a3);
        s24z = vertices(3,a2)-vertices(3,a4);
        h = max([s31x.^2+s31y.^2+s31z.^2;s21x.^2+s21y.^2+s21z.^2;...
            s41x.^2+s41y.^2+s41z.^2;s23x.^2+s23y.^2+s23z.^2;...
            s43x.^2+s43y.^2+s43z.^2;s24x.^2+s24y.^2+s24z.^2]);
        
        invjac(:,1,1) = dcdx;
        invjac(:,1,2) = dedx;
        invjac(:,1,3) = dtdx;
        
        invjac(:,2,1) = dcdy;
        invjac(:,2,2) = dedy;
        invjac(:,2,3) = dtdy;
        
        invjac(:,3,1) = dcdz;
        invjac(:,3,2) = dedz;
        invjac(:,3,3) = dtdz;
        
        detjac = abs(volumes);
        
end

return