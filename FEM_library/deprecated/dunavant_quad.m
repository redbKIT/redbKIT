function [quad_points, w] = dunavant_quad(degree)
%DUNAVANT Dunavant quadrature for a triangle.
%   [QUAD_POINTS,W]=DUNAVANT(DEGREE) computes the weights W and the quadrature nodes QUAD_POINTS
%   of the Dunavant quadrature on the unit triangle. DEGREE is the degree of the
%   complete polynomial exactly integrated.
%
%   [QUAD_POINTS,W]=DUNAVANT(DEGREE,AX,AY) the same for a generic triangle with vertices
%   (AX(I),AY(I)), I = 1,2,3.
%
%   THIS FUNCTION IS DEPRECATED. USE quadrature instead.

warning('dunavant_quad is deprecated. Use quadrature instead.')

[quad_points, w] = quadrature(2, degree);


return
