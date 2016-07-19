function [mu_samples, mu_greedy_index] = FullFactorial_ParameterSpace(P, mu_min, mu_max, mu_bar, mu_samples_Dimension, logflag)
%FULLFACTORIAL_PARAMETERSPACE generates a tensor product grid

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 6 || isempty(logflag)
      logflag = zeros(1,P);
end

if length(mu_samples_Dimension)==1 && P > 1
      mu_samples_Dimension_tmp = mu_samples_Dimension;
      
      mu_samples_Dimension     = ones(1,P) * ceil(mu_samples_Dimension_tmp^(1/P));
end

mu_greedy_index = [];
for i = 1 : P
      
      if mu_samples_Dimension(i) > 1
            mu_greedy_index = [mu_greedy_index i];
      end
      
end

output_args = '';
input_args  = '';
for i = 1 : length(mu_greedy_index)
      
      if logflag(i)
            eval(['ti',num2str(i),' =  logspace(log10(mu_min(mu_greedy_index(',num2str(i),'))), log10(mu_max(mu_greedy_index(',num2str(i),'))), mu_samples_Dimension(mu_greedy_index(',num2str(i),')));']);
      else
            eval(['ti',num2str(i),' =  linspace(mu_min(mu_greedy_index(',num2str(i),')), mu_max(mu_greedy_index(',num2str(i),')), mu_samples_Dimension(mu_greedy_index(',num2str(i),')));']);
      end
      
      input_args  = strcat(input_args,'ti',num2str(i),',');
      output_args = strcat(output_args,'q',num2str(i),',');
end
input_args  = input_args(1:end-1);
output_args = output_args(1:end-1);

eval(['[',output_args,'] = ndgrid(',input_args,');']);

for i = 1 : length(mu_greedy_index)
      eval(['q',num2str(i),' = q',num2str(i),'(:);']);
end

for i = 1:P
      mu_samples(:,i) = ones(length(q1),1)*mu_bar(i);
end

for i = 1 : length(mu_greedy_index)
      eval(['mu_samples(:,mu_greedy_index(',num2str(i),')) = q',num2str(i),';']);
end


return