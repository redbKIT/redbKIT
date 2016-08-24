function du0 = Initial_Velocity( x, y, z, t, param, d )
ang_vel = 300;

numNodes = length(x);


v1 = - ang_vel * sqrt(x.^2+y.^2) .* sin(atan2(y,x));

v2 = ang_vel * sqrt(x.^2+y.^2) .* cos(atan2(y,x));

v3 = 0.*x.*y.*z;



v = [v1'; v2'; v3'];

indx = find(isnan(v));
v(indx) = 0;


alpha = pi/3;
Rot_X_Matrix = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];

dU_n_tmp     = [v(1:numNodes)'; v(1+numNodes:2*numNodes)'; v(2*numNodes+1:end)'];

dU_n_tmp     = Rot_X_Matrix  * dU_n_tmp;
du0          = dU_n_tmp(d,:);

end