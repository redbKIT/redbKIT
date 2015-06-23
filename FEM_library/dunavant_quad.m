function [quad_points, w] = dunavant_quad(degree)
%DUNAVANT Dunavant quadrature for a triangle.
%   [QUAD_POINTS,W]=DUNAVANT(DEGREE) computes the weights W and the quadrature nodes QUAD_POINTS
%   of the Dunavant quadrature on the unit triangle. DEGREE is the degree of the
%   complete polynomial exactly integrated.
%
%   [QUAD_POINTS,W]=DUNAVANT(DEGREE,AX,AY) the same for a generic triangle with vertices
%   (AX(I),AY(I)), I = 1,2,3.
%

%  Author: F. Saleri 13-01-03.
%  Modified 27-02-07.

ax = [0 1 0];
ay = [0 0 1];
area = 0.5;


switch degree
    case 1
        w = 1;  x = sum(ax)/3; y = sum(ay)/3;
    case 2
        w = ones(1,3)/3;
        x = [2*ax(1)/3+ax(2)/6+ax(3)/6 ax(1)/6+2*ax(2)/3+ax(3)/6 ax(1)/6+ax(2)/6+2*ax(3)/3];
        y = [2*ay(1)/3+ay(2)/6+ay(3)/6 ay(1)/6+2*ay(2)/3+ay(3)/6 ay(1)/6+ay(2)/6+2*ay(3)/3];
    case 2.5
        % This is not a Dunavant quadrature
        w = ones(1,3)/3;
        x = [ax(1) ax(2) ax(3)];
        y = [ay(1) ay(2) ay(3)];
    case 3
        w = [-9/16 25/48 25/48 25/48];
        x = [sum(ax)/3 (3*ax(1)+ax(2)+ax(3))/5 (ax(1)+3*ax(2)+ax(3))/5 (ax(1)+ax(2)+3*ax(3))/5];
        y = [sum(ay)/3 (3*ay(1)+ay(2)+ay(3))/5 (ay(1)+3*ay(2)+ay(3))/5 (ay(1)+ay(2)+3*ay(3))/5];
        warning('The Dunavant formula of degree 3 presents negative weigths');
    case 4
        w = [0.223381589678011.*ones(1,3) 0.109951743655322.*ones(1,3)];
        p1 = 0.108103018168070; p2 = 0.445948490915965;
        x = [p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.816847572980459; p2 = 0.091576213509771;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
    case 5
        w = [9/40 0.132394152788506*ones(1,3) 0.125939180544827*ones(1,3)];
        x = sum(ax)/3;
        y = sum(ay)/3;
        p1 = 0.059715871789770;
        p2 = 0.470142064105115;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.797426985353087;
        p2 = 0.101286507323456;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
    case 6
        w = [0.116786275726379*ones(1,3) 0.050844906370207*ones(1,3) 0.082851075618374*ones(1,6)];
        p1 = 0.501426509658179;
        p2 = 0.249286745170910;
        x = [p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.873821971016996;
        p2 = 0.063089014491502;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.053145049844817;
        p2 = 0.310352451033784;
        p3 = 0.636502499121399;
        x = [x, p1*ax(1)+p2*ax(2)+p3*ax(3) p1*ax(1)+p3*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p3*ax(3), ...
            p2*ax(1)+p3*ax(2)+p1*ax(3) p3*ax(1)+p1*ax(2)+p2*ax(3) p3*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y,p1*ay(1)+p2*ay(2)+p3*ay(3) p1*ay(1)+p3*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p3*ay(3), ...
            p2*ay(1)+p3*ay(2)+p1*ay(3) p3*ay(1)+p1*ay(2)+p2*ay(3) p3*ay(1)+p2*ay(2)+p1*ay(3)];
    case 7
        w = [-0.149570044467682 0.175615257433208*ones(1,3) 0.053347235608838*ones(1,3) 0.077113760890257*ones(1,6)];
        x = sum(ax)/3;
        y = sum(ay)/3;
        p1 = 0.479308067841920;
        p2 = 0.260345966079040;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.869739794195568;
        p2 = 0.065130102902216;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.048690315425316;
        p2 = 0.312865496004874;
        p3 = 0.638444188569810;
        x = [x, p1*ax(1)+p2*ax(2)+p3*ax(3) p1*ax(1)+p3*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p3*ax(3), ...
            p2*ax(1)+p3*ax(2)+p1*ax(3) p3*ax(1)+p1*ax(2)+p2*ax(3) p3*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y,p1*ay(1)+p2*ay(2)+p3*ay(3) p1*ay(1)+p3*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p3*ay(3), ...
            p2*ay(1)+p3*ay(2)+p1*ay(3) p3*ay(1)+p1*ay(2)+p2*ay(3) p3*ay(1)+p2*ay(2)+p1*ay(3)];
    case 8
        w = [0.144315607677787 0.095091634267285*ones(1,3) 0.103217370534718*ones(1,3) 0.032458497623198*ones(1,3) ...
            0.027230314174435*ones(1,6)];
        x = sum(ax)/3;
        y = sum(ay)/3;
        p1 = 0.081414823414554;
        p2 = 0.459292588292723;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.658861384496480;
        p2 = 0.170569307751760;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.898905543365938;
        p2 = 0.050547228317031;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.008394777409958;
        p2 = 0.263112829634638;
        p3 = 0.728492392955404;
        x = [x, p1*ax(1)+p2*ax(2)+p3*ax(3) p1*ax(1)+p3*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p3*ax(3), ...
            p2*ax(1)+p3*ax(2)+p1*ax(3) p3*ax(1)+p1*ax(2)+p2*ax(3) p3*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y,p1*ay(1)+p2*ay(2)+p3*ay(3) p1*ay(1)+p3*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p3*ay(3), ...
            p2*ay(1)+p3*ay(2)+p1*ay(3) p3*ay(1)+p1*ay(2)+p2*ay(3) p3*ay(1)+p2*ay(2)+p1*ay(3)];
    case 9
        w = [0.097135796282799 0.031334700227139*ones(1,3) 0.077827541004774*ones(1,3) 0.079647738927210*ones(1,3) ...
            0.025577675658698*ones(1,3) 0.043283539377289*ones(1,6)];
        x = sum(ax)/3;
        y = sum(ay)/3;
        p1 = 0.020634961602525;
        p2 = 0.489682519198738;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.125820817014127;
        p2 = 0.437089591492937;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.623592928761935;
        p2 = 0.188203535619033;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.910540973211095;
        p2 = 0.044729513394453;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.036838412054736;
        p2 = 0.221962989160766;
        p3 = 0.741198598784498;
        x = [x, p1*ax(1)+p2*ax(2)+p3*ax(3) p1*ax(1)+p3*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p3*ax(3), ...
            p2*ax(1)+p3*ax(2)+p1*ax(3) p3*ax(1)+p1*ax(2)+p2*ax(3) p3*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y,p1*ay(1)+p2*ay(2)+p3*ay(3) p1*ay(1)+p3*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p3*ay(3), ...
            p2*ay(1)+p3*ay(2)+p1*ay(3) p3*ay(1)+p1*ay(2)+p2*ay(3) p3*ay(1)+p2*ay(2)+p1*ay(3)];
    case 10
        w = [0.090817990382754 0.036725957756467*ones(1,3) 0.045321059435528*ones(1,3) 0.072757916845420*ones(1,6) ...
            0.028327242531057*ones(1,6) 0.009421666963733*ones(1,6)];
        x = sum(ax)/3;
        y = sum(ay)/3;
        p1 = 0.028844733232685;
        p2 = 0.485577633383657;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.781036849029926;
        p2 = 0.109481575485037;
        x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.141707219414880;
        p2 = 0.307939838764121;
        p3 = 0.550352941820999;
        x = [x, p1*ax(1)+p2*ax(2)+p3*ax(3) p1*ax(1)+p3*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p3*ax(3), ...
            p2*ax(1)+p3*ax(2)+p1*ax(3) p3*ax(1)+p1*ax(2)+p2*ax(3) p3*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y,p1*ay(1)+p2*ay(2)+p3*ay(3) p1*ay(1)+p3*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p3*ay(3), ...
            p2*ay(1)+p3*ay(2)+p1*ay(3) p3*ay(1)+p1*ay(2)+p2*ay(3) p3*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.025003534762686;
        p2 = 0.246672560639903;
        p3 = 0.728323904597411;
        x = [x, p1*ax(1)+p2*ax(2)+p3*ax(3) p1*ax(1)+p3*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p3*ax(3), ...
            p2*ax(1)+p3*ax(2)+p1*ax(3) p3*ax(1)+p1*ax(2)+p2*ax(3) p3*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y,p1*ay(1)+p2*ay(2)+p3*ay(3) p1*ay(1)+p3*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p3*ay(3), ...
            p2*ay(1)+p3*ay(2)+p1*ay(3) p3*ay(1)+p1*ay(2)+p2*ay(3) p3*ay(1)+p2*ay(2)+p1*ay(3)];
        p1 = 0.009540815400299;
        p2 = 0.066803251012200;
        p3 = 0.923655933587500;
        x = [x, p1*ax(1)+p2*ax(2)+p3*ax(3) p1*ax(1)+p3*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p3*ax(3), ...
            p2*ax(1)+p3*ax(2)+p1*ax(3) p3*ax(1)+p1*ax(2)+p2*ax(3) p3*ax(1)+p2*ax(2)+p1*ax(3)];
        y = [y,p1*ay(1)+p2*ay(2)+p3*ay(3) p1*ay(1)+p3*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p3*ay(3), ...
            p2*ay(1)+p3*ay(2)+p1*ay(3) p3*ay(1)+p1*ay(2)+p2*ay(3) p3*ay(1)+p2*ay(2)+p1*ay(3)];
    otherwise
        error('Dunavant quadrature is limited to degree 10');
end

w = w*area;

quad_points = [x;y];

return