function [pmu] = Device_Heat_DefGeo(vertices, mu)

d      = mu(1);

a_0 = 1;
b_0 = 0;

a_1 = (1 + d  - 2/3)/(1- 2/3);
b_1 = 2/3*(1-a_1);


x = vertices(1,:);
y = vertices(2,:);

pmu(2,:)  = y;

pmu(1,:)  =   (a_0 * x + b_0) .* (x<=2/3) ...
            + (a_1 * x + b_1) .* (x>2/3);

end