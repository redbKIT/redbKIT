function [MESH] = FSI_InterfaceMap(DATA, MESH)
%FSI_INTERFACEMAP preprocessing function for FSI solver
%
%   [MESH] = FSI_INTERFACEMAP(DATA, MESH)
%   Generates mappings from solid to fluid interface dofs and viceversa.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

dim            = MESH.dim;
nodesS         = MESH.Solid.nodes;
boundariesS    = MESH.Solid.boundaries;
nodesF          = MESH.Fluid.nodes;
boundariesF    = MESH.Fluid.boundaries;
flag_interface = DATA.Fluid.flag_FSinterface;
flag_ALE       = DATA.Fluid.flag_ALE_fixed;
bcrow          = MESH.Fluid.bc_flag_row;
nbn            = MESH.Fluid.numBoundaryDof;
dirichletF     = MESH.Fluid.Dirichlet_dof_c;
dirichletS     = MESH.Solid.Dirichlet_dof_c;

for k = 1 : dim
    
    interfaceS_dofs{k} = [];
    for j = 1 : length(flag_interface{k})
        interfaceS_faces =  find(boundariesS(bcrow,:) == flag_interface{k}(j));
        tmp = boundariesS(1:nbn,interfaceS_faces);
        interfaceS_dofs{k}  = [interfaceS_dofs{k}; tmp(:)];
    end
    interfaceS_dofs{k}  = unique(interfaceS_dofs{k}(:));
    
    interfaceS_dofs{k}  = setdiff(interfaceS_dofs{k}, intersect(interfaceS_dofs{k},dirichletS{k}));
    
    
    interfaceF_dofs{k} = [];
    for j = 1 : length(flag_interface{k})
        interfaceF_faces =  find(boundariesF(bcrow,:) == flag_interface{k}(j));
        tmp = boundariesF(1:nbn,interfaceF_faces);
        interfaceF_dofs{k}  = [interfaceF_dofs{k}; tmp(:)];
    end
    interfaceF_dofs{k}  = unique(interfaceF_dofs{k}(:));
    
    interfaceF_dofs{k}  = setdiff(interfaceF_dofs{k}, intersect(interfaceF_dofs{k},dirichletF{k}));
    
    Interface_SFmap{k} = zeros(length(interfaceS_dofs{k}),1);
    Interface_FSmap{k} = zeros(length(interfaceF_dofs{k}),1);
   
    tmp = zeros(length(interfaceS_dofs{k}),1);
    parfor i = 1 : length(interfaceS_dofs{k})
        [~, iF] = utility(i, nodesS, interfaceS_dofs{k}, nodesF, interfaceF_dofs{k});
        
        tmp(i) = iF;
    end
    
    for i = 1 : length(interfaceS_dofs{k})        
        Interface_SFmap{k}(i)        =  tmp(i);
        
        Interface_FSmap{k}( tmp(i) ) = i;
    end
end

for k = 1 : dim
      % Computes the Dirichlet sides of the domain
      nDir = length(flag_ALE{k});
      Dirichlet_side_tmp = [];
      for j = 1 : nDir
          
            Dirichlet_side_j   = find(boundariesF(bcrow,:) == flag_ALE{k}(j));
            Dirichlet_side_tmp  = [Dirichlet_side_tmp, Dirichlet_side_j];
            
      end
      Dirichlet_side{k} = unique(Dirichlet_side_tmp);
      
      if ~isempty(Dirichlet_side{k})
            ALE_dirichlet{k} = boundariesF(1:nbn,Dirichlet_side{k});
            ALE_dirichlet{k} = unique(ALE_dirichlet{k}(:));
      else
            ALE_dirichlet{k} = [];
      end
end

MESH.Fluid.dof_interface =  interfaceF_dofs;
MESH.Solid.dof_interface =  interfaceS_dofs;
MESH.Interface_SFmap =  Interface_FSmap; % solid to fluid map
MESH.Interface_FSmap =  Interface_SFmap; % fluid to solid map
MESH.ALE_dirichlet   =  ALE_dirichlet;

return

function [i, iF] = utility(i, verticesS, interfaceS_dofs, verticesF, interfaceF_dofs)

iS_coord = verticesS(:,interfaceS_dofs(i));

iF = find( ismember(verticesF(:,interfaceF_dofs)', iS_coord', 'rows') );

return
