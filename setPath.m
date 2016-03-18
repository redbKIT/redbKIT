%SETPATH add the subfolder FEM_library and RB_library to the path 

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath(genpath(strcat(pwd,sslash,'RB_library')));
addpath(genpath(strcat(pwd,sslash,'FEM_library')));

fprintf('\n------------------------------------------')
fprintf('\n      *** Welcome to redbKIT v2.0 ***')
fprintf('\n------------------------------------------\n\n')
