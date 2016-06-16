function [phi, dphix, dphiy] = fem_basis2D(fem, x, y, varargin)
%FEM_BASIS2D finite element basis functions.
%       F. Saleri 13-01-03,  F. Negri 22.11.2014


warning('fem_basis2D is deprecated. Use fem_basis instead.')

[phi, dphi] = fem_basis(2, fem, [x;y], varargin{:});

dphix = dphi(:,:,1);
dphiy = dphi(:,:,2);

return
