function map = metis_to_matlab(AdjacencyMatrix, n_subdom, contiguous)
%METIS_TO_MATLAB wrapper to Metis multiway partition
%
%   map = METIS_TO_MATLAB(AdjacencyMatrix, n_subdom, contiguous)
%
%   Requires the adjacency matrix of the graph and the number of desired
%   subdomains. If contiguous = 1, metis returns a contiguous partition.
%
%   See also geometric_domain_decomposition
%
%   Inspired by metisdice (part of MESHPART toolbox) by John Gilbert
% 
%   Requires metis 5.0 and metismex interface. 
%   To compile metismex see the link below:
%   http://dgleich.wordpress.com/2012/05/22/want-to-use-metis-5-0-with-matlab-try-the-new-metismex/

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin == 2 || isempty(contiguous) || ~contiguous
    
      map = metismex('PartGraphKway', AdjacencyMatrix, n_subdom);
      
elseif nargin == 3
    
      wtflag.contig = 1;
      map = metismex('PartGraphKway', AdjacencyMatrix, n_subdom, wtflag);
      
end


