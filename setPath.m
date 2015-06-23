%SETPATH add the subfolder FEM_library and RB_library to the path 

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath(genpath(strcat(pwd,sslash,'RB_library')));
addpath(genpath(strcat(pwd,sslash,'FEM_library')));

fprintf('\n------------------------------------------')
fprintf('\n      *** Welcome to redbKIT ***')
fprintf('\n------------------------------------------\n\n')

%% Mex assembler
% To speed up the FEM assembly you should "mexify" the function 
% ADR_mex_assembler.m by running the following command in the Matlab prompt
% coder -build FEM_library/ADR_assembly.prj
% Then replace ADR_mex_assembler with ADR_mex_assembler_mex at line 109
% of Assembler_2D

%% To generate HTML documentation
% Download, extract and add to the path M2HTML:
%
% url_M2HTML = 'http://www.artefact.tk/software/matlab/m2html/m2html.zip';
% unzip(url_M2HTML);
% addpath(genpath(strcat(pwd,sslash,'m2html')));
%
% Then run the following command
% m2html('mfiles',{'RB_library', 'RB_library/RBF_interpolation' 'FEM_library'}, ...
%   'htmldir','Documentation', 'recursive','off', 'global','on',...
%   'template','frame', 'index','menu', 'graph','on');

