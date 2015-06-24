function [Pv] = Gram_Schmidt_orth(Z, v, X)
%GRAM_SCHMIDT_ORTH performs Gram-Schmidt orthonormalization process 
%
%   [PV] = GRAM_SCHMIDT_ORTH(Z, V) performs Gram-Schmidt orthonormalization
%   process on the vector V with respect to the basis Z and the euclidean
%   inner product. (Not that Z can be empty.)
%
%   [PV] = GRAM_SCHMIDT_ORTH(Z, V, X) performs Gram-Schmidt orthonormalization
%   process on the vector V with respect to the basis Z and the inner 
%   product X.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch> 

if nargin < 3 || isempty(X)
    X = speye(length(v), length(v));
end


if ~isempty(Z)
    
    alpha = Z'*(X*v);
    
    Pv    = v - Z*alpha;
    
    % cure for loss of orthogonality
    if norm(Pv) < 0.7*norm(v)
        alpha = Z'*(X*Pv);
        Pv    = Pv - Z*alpha;
    end
    
else
    
    Pv    = v;
    
end

if ~(norm(Pv) == 0)
    
    Pv = Pv/sqrt(dot(Pv,(X*Pv)));
    
end

end
