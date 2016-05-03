function DATA = CFD_read_DataFile(data_file, dim)
%CFD_READ_DATAFILE data_file parser
%
%   DATA = CFD_READ_DATAFILE(DATA_FILE) read the file specified by the string
%   DATA_FILE and put the fields values into the struct DATA

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 2 || isempty(dim)
    dim = 2;
end

%% Set Default values
for d = 1 : dim
    DATA.flag_dirichlet{d} = [];
    DATA.flag_neumann{d}   = [];
    DATA.flag_pressure{d}     = [];
end
        
switch dim
    
    case 2
        
        DATA.bcDir          = @(x,y,t,param)(0.*x.*y);
        DATA.bcNeu          = @(x,y,t,param)(0.*x.*y);
        DATA.force          = @(x,y,t,param)(0.*x.*y);
                            
    case 3
        
        DATA.bcDir          = @(x,y,z,t,param)(0.*x.*y);
        DATA.bcNeu          = @(x,y,z,t,param)(0.*x.*y);
        DATA.force          = @(x,y,z,t,param)(0.*x.*y);
end

DATA.kinematic_viscosity  = 1;
DATA.density              = 1;

%% Read data_file and put problem-data into the DATA struct

eval(data_file);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)
    
    eval(['DATA.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
    
end

[ DATA ] = dataParser( DATA );


end