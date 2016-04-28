function STR_export_solution( varargin )
%see CSM_EXPORT_SOLUTION

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

warning('STR_export_solution: deprecated function, please use CSM_export_solution instead');

CSM_export_solution(varargin{:});

end
