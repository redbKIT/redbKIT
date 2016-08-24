function [] = test( n_subdomains )

if nargin < 1
    n_subdomains = 4;
end

dim = 2;

%% build P1 mesh
[vertices, boundaries, elements] = initmesh('mesh_square','Jiggle','minimum','Hgrad',1.01,'Hmax',0.5);

%% refine mesh
for i = 1 : 5
    [vertices, boundaries, elements] = refinemesh('mesh_square',vertices, boundaries, elements);
end

fprintf('\n **** TEST METIS INSTALLATION ****\n');
fprintf(' * Spatial Dimension     = %d \n',dim);
fprintf(' * Number of Vertices    = %d \n',size(vertices,2));
fprintf(' * Number of Elements    = %d \n',size(elements,2));
fprintf(' * Number of Subdomains  = %d \n',n_subdomains);
fprintf('-------------------------------------------\n');

fprintf('\n Computing Adjacency Matrix (graph) of the mesh ...')
time_a = tic;
A   = compute_adjacency(vertices, elements, dim);
time_a = toc(time_a);
fprintf(' done in %f s', time_a);

fprintf('\n Partitioning the mesh ...')
time_p = tic;
map = metis_to_matlab(A, n_subdomains, 1);
time_p = toc(time_p);
fprintf(' done in %f s\n', time_p);

% To generate vtk of the decomposition uncomment the following two lines
% [~,~,~] = mkdir('Figures');
% geometric_domain_decomposition(vertices, elements, dim, n_subdomains, 1, 1, ['Figures/Subs',num2str(n_subdomains)], elements);


end