function [detjac, invjac, h] = geotrasf2D(vertices,elements)
%GEOTRASF2D linear two dimensional transformation

warning('geotrasf2D is deprecated. Use geotrasf instead.')

[detjac, invjac, h] = geotrasf(2, vertices, elements);

return
