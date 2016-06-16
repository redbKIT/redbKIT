function DATA = read_DataFile(data_file, dim, param)
%READ_DATAFILE data_file parser
%
%   DATA = READ_DATAFILE(DATA_FILE) read the file specified by the string
%   DATA_FILE and put the fields values into the struct DATA

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 2 || isempty(dim)
    dim = 2;
end

if nargin < 3 || isempty(param)
    param = [];
end

%% Set Default values
switch dim
    
    case 2
        DATA.flag_dirichlet = [];
        DATA.flag_neumann   = [];
        DATA.flag_robin     = [];
        
        DATA.bcDir          = @(x,y,t,param)(0.*x.*y);
        DATA.bcNeu          = @(x,y,t,param)(0.*x.*y);
        DATA.bcRob_alpha    = @(x,y,t,param)(0.*x.*y);
        DATA.bcRob_fun      = @(x,y,t,param)(0.*x.*y);
        
        DATA.force          = @(x,y,t,param)(0.*x.*y);
        
        DATA.diffusion      = @(x,y,t,param)(0.*x.*y);
        
        for i = 1 : dim
            DATA.transport{i}   = @(x,y,t,param)(0.*x.*y);
        end
        DATA.reaction       = @(x,y,t,param)(0.*x.*y);
        
        DATA.uexact         = @(x,y,t,param)(0.*x.*y);
        DATA.uxexact        = @(x,y,t,param)(0.*x.*y);
    
    
    case 3
        
        DATA.flag_dirichlet = [];
        DATA.flag_neumann   = [];
        DATA.flag_robin     = [];
        
        DATA.bcDir          = @(x,y,z,t,param)(0.*x.*y);
        DATA.bcNeu          = @(x,y,z,t,param)(0.*x.*y);
        DATA.bcRob_alpha    = @(x,y,z,t,param)(0.*x.*y);
        DATA.bcRob_fun      = @(x,y,z,t,param)(0.*x.*y);
        
        DATA.force          = @(x,y,z,t,param)(0.*x.*y);
        
        DATA.diffusion      = @(x,y,z,t,param)(0.*x.*y);
        
        for i = 1 : dim
            DATA.transport{i}   = @(x,y,z,t,param)(0.*x.*y);
        end
        DATA.reaction       = @(x,y,z,t,param)(0.*x.*y);
        
        DATA.uexact         = @(x,y,z,t,param)(0.*x.*y);
        DATA.uxexact        = @(x,y,z,t,param)(0.*x.*y);
        
end

%% Read data_file and put problem-data into the DATA struct

eval(data_file);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)
    
    eval(['DATA.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
    
end

[ DATA ] = dataParser( DATA );

end