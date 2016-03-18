function pp = RBF_perturbnodes_func(mu,vertices)
%RBF_PERTURBNODES_FUNC
%
%   def_vertices = RBF_PERTURBNODES_FUNC(mu,vertices)
%
%   Author: F. Negri (federico.negri@epfl.ch) 2013-2015
%   Copyright (C) Federico Negri, CMCS, EPFL




% Set shape parameters
for k = 1 : length(mu)
    eval(['mu' num2str(k) ' = mu(' num2str(k) ');']);
end
for k = length(mu)+1 : 8
    eval(['mu' num2str(k) ' = 0;']);
end


% Set nodes
x = vertices(1,:);
y = vertices(2,:);


pp(1,:) = T_1(mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8, x, y);
pp(2,:) = T_2(mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8, x, y);


return