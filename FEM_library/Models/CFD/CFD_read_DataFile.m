function DATA = CFD_read_DataFile(data_file, dim, param)
%CFD_READ_DATAFILE data_file parser
%
%   DATA = CFD_READ_DATAFILE(DATA_FILE, DIM, PARAM) read the file specified by the string
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
    DATA.flag_pressure{d}    = [];
    DATA.flag_ring{d}        = [];
end
DATA.flag_resistance  = [];
DATA.flag_absorbing   = [];

switch dim
    
    case 2
        
        for d = 1 : dim
            DATA.bcDir{d}          = @(x,y,t,param)(0.*x.*y);
            DATA.bcNeu{d}          = @(x,y,t,param)(0.*x.*y);
            DATA.force{d}          = @(x,y,t,param)(0.*x.*y);
        end
        
    case 3
        
        for d = 1 : dim
            DATA.bcDir{d}          = @(x,y,z,t,param)(0.*x.*y);
            DATA.bcNeu{d}          = @(x,y,z,t,param)(0.*x.*y);
            DATA.force{d}          = @(x,y,z,t,param)(0.*x.*y);
        end
end

DATA.dynamic_viscosity    = 1;
DATA.density              = 1;

%% Read data_file and put problem-data into the DATA struct

eval(data_file);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)
    
    eval(['DATA.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
    
end

[ DATA ] = dataParser( DATA );

end