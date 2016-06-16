function RBF_data = RBF_setup(x, y, RBF_function_name, constant)
%RBF_SETUP computes RBF coefficient 
%
%   RBF_DATA = RBF_SETUP(X, Y, RBF_FUNCTION) given a matrix X of size dimX
%   x num_interpolation_points, a row vector Y of size 1 x
%   num_interpolation_points, and a string RBF_FUNCTION specifying the
%   type of RBF to be used, returns a struct containing the RBF
%   interpolation coefficients

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


if nargin < 3 || isempty(RBF_function_name)
    RBF_function_name = 'gaussian';
end
    

[dimX, num_int_p] = size(x);

if (length(y)~=num_int_p)
  error('x and y should have the same number of columns');
end

if (size(y,1)~=1)
  error('y should be a row vector');
end

RBF_data.x    = x;
RBF_data.y    = y;

RBF_data.RBF_function_type  = RBF_function_name;

if nargin < 4 || isempty(constant)
    RBF_data.constant      = (prod(max(x,[],2)-min(x,[],2))/num_int_p)^(1/dimX); %approx. average distance between the nodes
else
    RBF_data.constant = constant;
end

switch RBF_function_name
    
      case 'cubic'
        RBF_data.RBF_function   = @(r,c)(r.*r.*r);
        
      case 'multiquadric'
        RBF_data.RBF_function   = @(r,c)sqrt(1+r.*r/(c*c));
        
      case 'thinplate'
        RBF_data.RBF_function   = @(r,c)r.*r.*log(r+1);
        
      case 'gaussian'
        RBF_data.RBF_function   = @(r,c)exp(-0.5*r.*r/(c*c));
        
    otherwise
        warning('RBF_function set to gaussian')
        RBF_data.RBF_function   = @(r,c)exp(-0.5*r.*r/(c*c));
end

A       = zeros(num_int_p,num_int_p);
for i = 1 : num_int_p
    for j = 1 : i
        r      = norm(x(:,i)-x(:,j));
        A(i,j) = RBF_data.RBF_function(r, RBF_data.constant);
        A(j,i) = A(i,j);
    end
end


P   = [ones(num_int_p,1) x'];
A   = [ A                P
        P' zeros(dimX+1,dimX+1)];
  
b   = [y'; zeros(dimX+1, 1)];                       

RBF_data.coeff = A \ b;

end

