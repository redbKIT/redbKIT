function A = GlobalAssemble( i, j, s, m, n )
%GLOBALASSEMBLE Create sparse matrix from elemental contributions
% 
%   A = GLOBALASSEMBLE(i,j,s,m,n) uses vectors i, j, and s to generate an
%   m-by-n sparse matrix such that A(i(k),j(k)) = s(k). Vectors i, j, and s 
%   are all the same length.  Any elements of s that are zero are ignored, 
%   along with the corresponding values of i and j.  Any elements of s that 
%   have duplicate values of i and j are added together.  The argument s and one of the
%   arguments i or j may be scalars, in which case the scalars are expanded
%   so that the first three arguments all have the same length.
%
%   If the FAST package is available, then GLOBALASSEMBLE calls fsparse,
%   otherwise it calls Matlab built-in sparse function.
% 
%   The sintax of GLOBALASSEMBLE is the same as that of sparse and fsparse;
%   please type help sparse or help sparse. For instance
%
%   The FAST package is available at 
%       http://user.it.uu.se/~stefane/freeware
%
%   See also the companion paper
%   S. Engblom, D. Lukarski: Fast Matlab compatible sparse assembly on multicore computers, 
%   in Parallel Comput. 56:1--17 (2016)
%
%   See also sparse, fsparse.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if exist('fsparse', 'file') == 3
    
    A = fsparse( i, j, s, [m, n] );
    
else
    
    A = sparse( i, j, s, m, n );
    
end

end