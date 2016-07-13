function [] = test(fem, datafile)

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

dim = 3;
h   = [2/4 2/9];

%% solve for different level of refinement
for i = 1 : 2
    
    %% load P1 mesh
    [vertices, boundaries, elements] = msh_to_Mmesh(strcat('gmsh/Cubeh',num2str(i)), dim);
    
    %% Solve   
    [~, ~, ~, ~, errorL2(i), errorH1(i)]  = Elliptic_Solver(dim, elements, vertices, boundaries, fem, datafile);
        
end

%% Plot convergence order
r_FE = str2double(fem(2));

figure
loglog(h,errorL2,'-ob','LineWidth',2)
hold on
loglog(h,errorH1,'-or','LineWidth',2)
loglog(h, errorH1(1)/2/h(1)*h.^r_FE, 'r--', h, errorL2(1)/2/(h(1)^2)*h.^(r_FE+1), 'b--')
legend(['errL2',fem],['errH1',fem], ['h^',num2str(r_FE)],  ['h^',num2str(r_FE+1)])
title(['Convergence ',fem,' Finite elements'])
xlim([min(h) max(h)])
grid on

end