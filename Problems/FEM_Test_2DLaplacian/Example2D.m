%EXAMPLE2D shows how to solve a problem with the FEM library. Also performs
%convergence analysis.

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

clc
clear all

fem      =  'P1'; 

%% build P1 mesh
[vertices, boundaries, elements] = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax',0.5);


%% solve for different level of refinement
for i = 1 : 5
    
    %% refine mesh
    [vertices, boundaries, elements] = refinemesh('mesh_square',vertices, boundaries, elements);
    h(i) = 1/(2^(i+1));
    
    %% Solve   
    [U, FE_SPACE, MESH, DATA, errorL2(i), errorH1(i)]  = Elliptic2D_Solver(elements, vertices, boundaries, fem, 'Dirichlet_data');
    %[U, FE_SPACE, MESH, DATA, errorL2(i), errorH1(i)] = Elliptic2D_Solver(elements, vertices, boundaries, fem, 'DirichletNeumann_data');
    %[U, FE_SPACE, MESH, DATA, errorL2(i), errorH1(i)] = Elliptic2D_Solver(elements, vertices, boundaries, fem, 'DirichletRobin_data');
    
    %% Visualize solution
    figure
    pdeplot(MESH.vertices,[],MESH.elements(1:3,:),'xydata',U(1:MESH.numVertices),'xystyle','interp',...
       'zdata',U(1:MESH.numVertices),'zstyle','continuous',...
       'colorbar','on','mesh','on');
    colormap(jet);
    lighting phong

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

