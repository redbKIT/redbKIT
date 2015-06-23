function [f] = nonaffine_source(x, mu)

f =    exp(-((x(:,1)-mu(1)).^2 + (x(:,2)-mu(2)).^2)/(0.25^2) );

end