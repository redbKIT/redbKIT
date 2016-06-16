%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

mu_test_g = [ 0.0 0.0 0.0 0.0;...
              0.03 0.01 0.02 0.03;...
             -0.02 0.01 0.015 0.03;...
             -0.03 -0.02 -0.03 -0.028];

GAMMA_IN_side = find(ROM.MESH.boundaries(5,:)==2);

GAMMA_IN_vert = ROM.MESH.boundaries(1:ROM.MESH.numBoundaryDof, GAMMA_IN_side);

GAMMA_IN_vert = unique(GAMMA_IN_vert(:)); 

handles = figure;
for jj = 1 : 4
      
      mu_test_Dimension = [150 1 1 1 1];
      [mu_test]         = FullFactorial_ParameterSpace(FOM.P, FOM.mu_min, FOM.mu_max, [500 mu_test_g(jj,:)], mu_test_Dimension);
            
      def_vertices = RBF_MeshDeformation(mu_test_g(jj,:), ROM.MESH.vertices);
      
      M_in = AssembleMass1D( FOM.FE_SPACE.fem, ROM.MESH.boundaries, def_vertices, ROM.MESH.numNodes, GAMMA_IN_side, GAMMA_IN_vert);
      
      Jrom     = zeros(size(mu_test,1),1);
      Jfom     = zeros(size(mu_test,1),1);

      for i = 1 : size(mu_test,1)
          
          Jrom(i)       = output_IRI(mu_test(i,:), ROM, GAMMA_IN_vert, M_in, ROM.N);
          
          [Uh]          = Elliptic_Solver(2, ROM.MESH.elements, def_vertices, ROM.MESH.boundaries, FOM.FE_SPACE.fem, 'horn_data', mu_test(i,:));
          
          J             = 1/0.05 * ones(1,length(GAMMA_IN_vert))*M_in*Uh(GAMMA_IN_vert);
          
          Jfom(i) = abs(J-1);
          
      end
      
      % Visualize
      freq = mu_test(:,1);
      
      subplot(2,2,jj)
      plot(freq, Jrom,'-b','LineWidth',3)
      hold on
      plot(freq, Jfom,'-r','LineWidth',3)
      legend('ROM','FOM')
      xlabel('Frequency (mu1) [Hz]')
      ylabel('Reflection spectrum')
      grid on
      xlim([FOM.mu_min(1) FOM.mu_max(1)])
      ylim([0 1])
      
      title_str = sprintf('%2.2f ',mu_test(1,2:end));
      title_str = strcat('mu_g = [',title_str,']');
      
      title(title_str);%['mu = ',num2str(mu_test(1,:))])
      
end

set(findall(handles,'-property','FontSize'),'FontSize',12)
saveas(handles,'Figures/comparison_IRI_FOM_ROM','fig');

