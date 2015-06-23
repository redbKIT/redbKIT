function [detjac, invjac, h] = geotrasf2D(vertices,elements)
%GEOTRASF2D linear two dimensional transformation
%
%   F. Saleri 24-08-01, F. Negri 18.11.2014

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

return