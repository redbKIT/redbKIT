function [ R ] = FSI_overlapping_DD( MESH, n_subdom, overlap, n_aggregates )
%FSI_OVERLAPPING_DD builds restriction operators associated to mesh decompositions
%for FSI problems in the two-fields condensed formulation
%
%   [ R ] = FSI_OVERLAPPING_DD( MESH , N_SUBDOM, OVERLAP )
%   given a MESH data structure (see also buildMESH.m and FSIt_Solver.m), the number of
%   subdomains N_SUBDOM and the overlap level OVERLAP (>= 1), returns a cell
%   array R of length N_SUBDOM. Each element of R (say R{i}) contains a
%   vector with the indices of the DOFs pertaining to the i-th subdomain.
%
%   see also geometric_domain_decomposition, metis_to_matlab,
%   AS_Preconditioner, FSIt_Solver

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


if nargin < 4 || isempty(n_aggregates)
    compute_coarse_aggregates = false;
else
    compute_coarse_aggregates = true;
end

dim      = MESH.dim;
nln      = MESH.Fluid.numElemDof;

F_vertices = MESH.Fluid.vertices ;

F_elements = MESH.Fluid.elements(1:dim+1, :);% P1 elements

F_vertices_fem = MESH.Fluid.nodes;

F_elements_fem = MESH.Fluid.elements(1:nln, :);

F_internal_vertices = MESH.Fluid.internal_dof_c;

F_II = MESH.Fluid.internal_dof(MESH.Fluid.II);

S_vertices = MESH.Solid.vertices;
S_elements = MESH.Solid.elements(1:dim+1, :);% P1 elements

S_vertices_fem = MESH.Solid.nodes;
S_elements_fem = MESH.Solid.elements(1:nln, :);

S_internal_vertices =  MESH.Solid.internal_dof_c;
S_II = MESH.Solid.internal_dof(MESH.Solid.II);

F_dof_interface = MESH.Fluid.dof_interface{1};
S_dof_interface = MESH.Solid.dof_interface{1};
F_Gamma = MESH.Fluid.internal_dof(MESH.Fluid.Gamma);
Interface_FSmap = MESH.Interface_FSmap{1};


%% Create Fluid+Solid+Interface mesh data structure
nov_S                = size(S_vertices_fem,2);
nov_F                = size(F_vertices_fem,2);

nov_SP1              = size(S_vertices,2);
nov_FP1              = size(F_vertices,2);

S_not_interface_dofs   = setdiff([1:nov_S],   S_dof_interface) ;
S_not_interface_dofsP1 = setdiff([1:nov_SP1], S_dof_interface) ;

S_dof_interfaceP1        = intersect(S_dof_interface, [1:nov_SP1]);
F_dof_interfaceP1        = intersect(F_dof_interface, [1:nov_FP1]);

F_elements_fem = F_elements_fem(1:nln,:);
S_elements_fem = S_elements_fem(1:nln,:);

F_elements = F_elements(1:dim+1,:);
S_elements = S_elements(1:dim+1,:);

shift_index                        = 0*[1:nov_S];
shift_index(S_not_interface_dofs)  = nov_F+[1:length(S_not_interface_dofs)];
shift_index(S_dof_interface)       = F_dof_interface(Interface_FSmap);

shift_indexP1                          = 0*[1:nov_SP1];
shift_indexP1(S_not_interface_dofsP1)  = nov_FP1+[1:length(S_not_interface_dofsP1)];
shift_indexP1(S_dof_interfaceP1)       = F_dof_interfaceP1(Interface_FSmap(1:length(F_dof_interfaceP1)));

vertices_fem       = [F_vertices_fem(1:dim, :)   S_vertices_fem(1:dim,S_not_interface_dofs)]; 
elements_fem       = [F_elements_fem(1:nln,:)    shift_index(S_elements_fem(1:nln, :))];
vertices           = [F_vertices_fem(1:dim, 1:nov_FP1)  S_vertices_fem(1:dim,S_not_interface_dofsP1)]; 
elements           = [F_elements_fem(1:dim+1,:)  shift_indexP1(S_elements_fem(1:dim+1, :))];

shift_back       = zeros(1, size(vertices_fem,2));
shift_back(nov_F+[1:length(S_not_interface_dofs)]) = S_not_interface_dofs;       
shift_back(F_dof_interface(Interface_FSmap)) = S_dof_interface;   

%% Partition the Fluid+Solid+Interface mesh
[subdom, ~, ~, A_elemC] = geometric_domain_decomposition(vertices, elements, dim, n_subdom, overlap, 1, 'Figures', elements_fem);


%% Translate the partition into the fluid and solid partitions, separately
for i = 1 : n_subdom
   
    subdomF{i} = intersect(1:nov_F, subdom{i});
    
    subdomS{i} = shift_back ( intersect(union(F_dof_interface, (1+nov_F):size(vertices_fem,2)), subdom{i}) );

end

% for i = 1 : n_subdom
%     subdomF2{i} =  zeros(size(nov_F,2),1);
%     subdomF2{i}(subdomF{i}) = subdomF{i};
%     
%     subdomS2{i} =  zeros(size(nov_S,2),1);
%     subdomS2{i}(subdomS{i}) = subdomS{i};
% end

%% restrict subdomains to internal vertices for each component of velocity and pressure
F_component = (dim + 1);
for k = 1 : F_component
      for i = 1 : n_subdom
            subdom_I{i,k} = intersect(subdomF{i}, F_internal_vertices{k});%nonzeros(subdomF2{i}(F_internal_vertices{k}));%
            if size(subdom_I{i,k},1) > size(subdom_I{i,k},2)
                  subdom_I{i,k} = subdom_I{i,k}';
            end
      end
end

S_component = dim;
for k = F_component+1 : S_component+F_component
      for i = 1 : n_subdom
              %subdom_I{i,k} = intersect(subdomS{i}, S_internal_vertices{k-F_component});
              subdom_I{i,k} = intersect(subdomS{i}, setdiff(S_internal_vertices{k-F_component},S_dof_interface));%nonzeros(subdomS2{i}(S_internal_vertices{k-F_component}));%
            if size(subdom_I{i,k},1) > size(subdom_I{i,k},2)
                  subdom_I{i,k} = subdom_I{i,k}';
            end
      end
end

%% build subdomains restriction/prolongation operators
nov_v    = size(F_vertices_fem,2);
nov_p    = size(F_vertices,2);
nov_totF = dim*nov_v + nov_p;
nov_totS = dim*nov_S;

nov_tot  = nov_totF+nov_totS;

I_all    = [F_II; F_Gamma; nov_totF+S_II];

for i = 1 : n_subdom
   
    subdom_I_all = [];
  
    for k = 1 : F_component
        subdom_I_all = [subdom_I_all nov_v*(k-1)+subdom_I{i,k}];
    end
    
    for k = F_component+1 : F_component+S_component
        subdom_I_all = [subdom_I_all nov_totF+nov_S*(k-F_component-1)+subdom_I{i,k}];
    end

    tmp =  zeros(nov_tot,1);
    tmp(subdom_I_all) = 1;
    tmp2  =  tmp(I_all);    
    R{i}  = find(tmp2);
    
end

check_n_subd = length(R);
fprintf('\n%d subdomains and restriction/prolongation operators built ---\n',check_n_subd);
        

%% build coarse aggregation restriction/prolongation operator

% if compute_coarse_aggregates
%     
%     [subdom_Coarse] = geometric_aggregates(A_elemC, vertices, elements, dim, n_aggregates, elements_fem);
%     
%     % Translate the partition into the fluid and solid partitions, separately
%     for i = 1 : n_aggregates
%         
%         subdom_CF{i} = intersect(1:nov_F, subdom_Coarse{i});
%         subdom_CS{i} = shift_back ( intersect(union(F_dof_interface, (1+nov_F):size(vertices_fem,2)), subdom_Coarse{i}) );
%         
%     end
%     
%     R{n_subdom+1}  = sparse((F_component+S_component)*n_aggregates, nov_tot);
%     
%     for i = 1 : n_aggregates
%         
%         for k = 1 : F_component
%             R{n_subdom+1}(i+n_aggregates*(k-1), nov_v*(k-1)+subdom_Coarse{i}) = 1;
%             %R{n_subdom+1}(k+i*(n_aggregates-1), nov_v*(k-1)+subdom_Coarse{i}) = 1;
%         end
%         for k = F_component+1 : F_component+S_component
%             R{n_subdom+1}(i+n_aggregates*(k-1), nov_totF+nov_S*(k-F_component-1)+subdom_Coarse{i}) = 1;
%             %R{n_subdom+1}(k+i*(n_aggregates-1), nov_totF+nov_S*(k-F_component-1)+subdom_Coarse{i}) = 1;
%         end
%         
%     end
%     
%     R{n_subdom+1} = R{n_subdom+1}(:,I_all);
%     
%     fprintf('\n%d Coarse aggregates computed ---\n', n_aggregates);
%     
% end

if compute_coarse_aggregates
    
    [subdom_CoarseF] = geometric_aggregates([], F_vertices, F_elements, dim, n_aggregates(1), F_elements_fem, 'FluidAggregates');
    [subdom_CoarseS] = geometric_aggregates([], S_vertices, S_elements, dim, n_aggregates(2), S_elements_fem, 'SolidAggregates');
    
    % Translate the partition into the fluid and solid partitions, separately
    for i = 1 : n_aggregates(1)
        subdom_CF{i} = intersect(1:nov_F, subdom_CoarseF{i});
    end
    
    for i = 1 : n_aggregates(2)
        subdom_CS{i} = intersect(1:nov_S, subdom_CoarseS{i});
    end
    
    R{n_subdom+1}  = sparse( F_component*n_aggregates(1)+S_component*n_aggregates(2), nov_tot);
    
    %     for i = 1 : n_aggregates
    %
    %         for k = 1 : F_component
    %             R{n_subdom+1}(i+n_aggregates*(k-1), nov_v*(k-1)+subdom_CF{i}) = 1;
    %             %R{n_subdom+1}(k+i*(n_aggregates-1), nov_v*(k-1)+subdom_Coarse{i}) = 1;
    %         end
    %         for k = F_component+1 : F_component+S_component
    %             R{n_subdom+1}(i+n_aggregates*(k-1), nov_totF+nov_S*(k-F_component-1)+subdom_CS{i}) = 1;
    %             %R{n_subdom+1}(k+i*(n_aggregates-1), nov_totF+nov_S*(k-F_component-1)+subdom_Coarse{i}) = 1;
    %         end
    %
    %     end
    
    for i = 1 : n_aggregates(1)
        
        for k = 1 : F_component
            R{n_subdom+1}(i+n_aggregates(1)*(k-1), nov_v*(k-1)+subdom_CF{i}) = 1;
            %R{n_subdom+1}(k+i*(n_aggregates-1), nov_v*(k-1)+subdom_Coarse{i}) = 1;
        end
        
    end
    
    for i = 1 : n_aggregates(2)
        
        for k = F_component+1 : F_component+S_component
            R{n_subdom+1}(i+n_aggregates(2)*(k-1-F_component)+n_aggregates(1)*F_component, nov_totF+nov_S*(k-F_component-1)+subdom_CS{i}) = 1;
            %R{n_subdom+1}(k+i*(n_aggregates-1), nov_totF+nov_S*(k-F_component-1)+subdom_Coarse{i}) = 1;
        end
        
    end
    
    R{n_subdom+1} = R{n_subdom+1}(:,I_all);
    
    fprintf('\n%d Coarse aggregates computed ---\n', n_aggregates);
    
end



return