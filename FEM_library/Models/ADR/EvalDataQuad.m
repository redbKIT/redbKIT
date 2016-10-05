function [f_quad] = EvalDataQuad(MESH, param, FE_SPACE, t, X_quad, f, index_subd)
%EVALDATAQUAD Evaluation of coefficients in the quadrature nodes

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch> 


switch class(f)
    
    case {'function_handle','inline'}
        
        switch MESH.dim
            case 2
                f_quad  = f(X_quad{1},X_quad{2},t,param);
            case 3
                f_quad  = f(X_quad{1},X_quad{2},X_quad{3}, t, param);
        end
        
    case 'double'
        
        if size(f,1)>1 && size(f,2)>1
            f_quad = f(index_subd,:);
            
        elseif length(f)==MESH.numNodes
            
            f_quad = zeros(MESH.numElem,FE_SPACE.numQuadNodes);
            for j = 1 : FE_SPACE.numElemDof
                i = MESH.elements(j,:);
                f_quad = f_quad + f(i)*FE_SPACE.phi(j,:);
            end
        end
end

end