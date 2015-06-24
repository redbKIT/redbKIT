function DATA = read_DataFile(data_file)
%READ_DATAFILE data_file parser
%
%   DATA = READ_DATAFILE(DATA_FILE) read the file specified by the string 
%   DATA_FILE and put the fields values into the struct DATA

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

%% Set Default values
DATA.flag_dirichlet = [];
DATA.flag_neumann   = [];
DATA.flag_robin     = [];

DATA.bcDir          = @(x,y,t,param)(0.*x.*y);
DATA.bcNeu          = @(x,y,t,param)(0.*x.*y);
DATA.bcRob_alpha    = @(x,y,t,param)(0.*x.*y);
DATA.bcRob_fun      = @(x,y,t,param)(0.*x.*y);

DATA.force          = @(x,y,t,param)(0.*x.*y);

DATA.diffusion      = @(x,y,t,param)(0.*x.*y);
DATA.transport{1}   = @(x,y,t,param)(0.*x.*y);
DATA.transport{2}   = @(x,y,t,param)(0.*x.*y);
DATA.reaction       = @(x,y,t,param)(0.*x.*y);

DATA.uexact         = @(x,y,t,param)(0.*x.*y);
DATA.uxexact        = @(x,y,t,param)(0.*x.*y);

%% Read data_file and put problem-data into the DATA struct

eval(data_file);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)
   
      eval(['DATA.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
      
end

end