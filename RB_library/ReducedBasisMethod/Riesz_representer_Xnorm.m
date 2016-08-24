function [R] = Riesz_representer_Xnorm(x, L, LT, p)
%RIESZ_REPRESENTER_XNORM given a vector x compute its "Riesz representer"
%R = X \ x.  
%
%   If a Cholesky factorization of X = L*LT (with reordering p) is provided,
%   then R = LT \ L \ x(p).
%
%   x can be either a vector or a "matrix of snapshots". 

%   Author: F. Negri (federico.negri@epfl.ch) 2013-2014
%   Copyright (C) Federico Negri, CMCS, EPFL

if nargin < 4 || isempty(p)
      p     = 1 : size(x,1);
end



if nargin == 2
      R = L \ x(p,:);
      
elseif nargin > 2
      R = LT \ ( L \ x(p,:));
end


if nargin == 4
      invp             = 0*p ;
      invp(p)          = 1:length(p);
      
      R = R(invp,:);
end

end
