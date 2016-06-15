%% Perform Error analysis
Ntest_sample  = 30;
[ Error ]     = error_analysis(FOM, ROM, Ntest_sample);


%% Compute Output of interest
GAMMA_IN_side = find(ROM.MESH.boundaries(5,:)==2);

GAMMA_IN_vert = ROM.MESH.boundaries(1:ROM.MESH.numBoundaryDof, GAMMA_IN_side);

GAMMA_IN_vert = unique(GAMMA_IN_vert(:)); 

M_in = AssembleMass1D( FOM.FE_SPACE.fem, ROM.MESH.boundaries, ROM.MESH.vertices, ROM.MESH.numNodes, GAMMA_IN_side, GAMMA_IN_vert);

mu_test_Dimension = 250;
mu_test           = linspace(FOM.mu_min, FOM.mu_max, mu_test_Dimension)';

Jrom     = zeros(size(mu_test,1),3);
Jfom     = zeros(size(mu_test,1),1);
rom_time = zeros(size(mu_test,1),3);
N_vec    = [20 30 ROM.N];


for j = 1 : 3
      parfor i = 1 : size(mu_test,1)
            tmp_time    = tic;
            
            Jrom(i,j)     = output_IRI(mu_test(i,:), ROM, GAMMA_IN_vert, M_in, N_vec(j));
            
            rom_time(i,j) = toc(tmp_time);
      end
end

parfor i = 1 : size(mu_test,1)
      
      tmp_time      = tic;
      
      Uh            = solve_HFsystem(FOM, mu_test(i,:));
      J             = 1/0.05 * ones(1,length(GAMMA_IN_vert))*M_in*Uh(GAMMA_IN_vert);

      Jfom(i) = abs(J-1);
      fom_time(i) = toc(tmp_time);
end


%% Visualize
freq = mu_test(:,1);

handle1=figure;
semilogy(freq, Jfom,'-r','LineWidth',3)
hold on
semilogy(freq, Jrom(:,1),'-g','LineWidth',2)
semilogy(freq, Jrom(:,2),'-k','LineWidth',2)
semilogy(freq, Jrom(:,3),'--b','LineWidth',2)
legend('IRI FOM','IRI ROM N=20','IRI ROM N=30',['IRI ROM N=',num2str(ROM.N)])
xlabel('Frequency \mu_1 [Hz]')
ylabel('Reflection spectrum')
grid on
xlim([FOM.mu_min(1) FOM.mu_max(1)])
set(findall(handle1,'-property','FontSize'),'FontSize',14)

saveas(handle1,'Figures/Reflection_spectrum_Comparison','fig');


fprintf('\n--------------------------------------------------------')
fprintf('\n Average ROM time = %2.2e s (N=20), %2.2e s (N=30), %2.2e s (N=50)',...
      mean(rom_time(:,1)),mean(rom_time(:,2)),mean(rom_time(:,3)))
fprintf('\n Average FOM time = %2.2e s',mean(fom_time))
fprintf('\n--------------------------------------------------------\n')