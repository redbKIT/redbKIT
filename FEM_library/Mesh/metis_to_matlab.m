function map = metis_to_matlab(AdjacencyMatrix, n_subdom, contiguous)
%MLIFE_METISDICE : Metis multiway partition
%
% map = metis_to_MLife(AdjacencyMatrix, n_subdom, contiguous)
%
% Requires the adjacency matrix of the graph and the number of desired
% subdomains. If contiguous = 1, metis returns a contiguous partition.
%
%
% Inspired by metisdice (part of MESHPART toolbox) by John Gilbert
% 
% Requires metis 5.0 and metismex interface. 
% To compile metismex see the link below:
% http://dgleich.wordpress.com/2012/05/22/want-to-use-metis-5-0-with-matlab-try-the-new-metismex/
%

%   Author: F. Negri (federico.negri@epfl.ch) 2013-2014
%   Copyright (C) Federico Negri, CMCS, EPFL
%

if nargin == 2 || isempty(contiguous) || ~contiguous
    
      map = metismex('PartGraphKway',AdjacencyMatrix,n_subdom);
      
elseif nargin == 3
    
      wtflag.contig = 1;
      map = metismex('PartGraphKway',AdjacencyMatrix,n_subdom,wtflag);
      
end


