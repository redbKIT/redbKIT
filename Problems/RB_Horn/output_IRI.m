function J = output_IRI(mu, ROM, GAMMA_IN_vert, M_gamma, N)
%OUTPUT_IRI computes the index of reflection intensity 

      Ur  = solve_RBsystem(ROM, mu, N);
      
      Ur_full_gamma = ROM.V(GAMMA_IN_vert,1:N)*Ur;
      
      a   = ones(1,length(GAMMA_IN_vert))*M_gamma*ones(1,length(GAMMA_IN_vert))';
      
      J   = 1/a * ones(1,length(GAMMA_IN_vert)) * (M_gamma*Ur_full_gamma);
      
      J   = abs(J-1);
      
end