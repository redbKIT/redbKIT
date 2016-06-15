function install_FAST(hasOpenMP)
%INSTALL_FAST installs FAST for sparse matrix assembly
%
%   INSTALL_FAST(hasOpenMP) if hasOpenMP = 1 enables OpenMP parallelism.
%   Default is hasOpenMP = 0.
%
%   For further details see:
%   http://user.it.uu.se/~stefane/freeware

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


FAST_license = sprintf('%s\n%s\n%s\n%s\n%s',...
    'You may download all software on this page and use, modify and redistribute it in any way you like.', ...
    'A redistributor must fully attribute the authorship and make a good effort to cite the original location of the software.', ...
    'A researcher making critical use of the software in research is requested to acknowledge this in publications related to the research.', ...
    'A company may use the code in software products provided that the original location and the author is clearly cited.', ...
    'All code provided here comes with absolutely no warranty and no support whatsoever is given');

fprintf('\n******************* FAST License by S. Engblom ********************');
fprintf('\n%s', FAST_license);
fprintf('\n\nFor more information see http://user.it.uu.se/~stefane/freeware');
fprintf('\n********************************************************************\n');

reply = input('\nDo you agree with the above license? (type Y or N):','s');

if strcmp(reply,'N')
    return;
end

if nargin < 1 || isempty(hasOpenMP)
    hasOpenMP = 0;
end

%% Download FAST
url_FAST = 'http://user.it.uu.se/cgi-bin/cgiwrap/~stefane/countFast.cgi?Fsparse_tests.tar';

fprintf('Download and untar parallel FAST from\n%s\n...\n',url_FAST);

untar(url_FAST);

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath(genpath(strcat(pwd,sslash,'Fast')));

cd Fast/

if hasOpenMP
    
    fprintf('\nCompiling FAST with openmp enabled\n');
    make('openmp',true,'fsparseonly',1);

else
    
    fprintf('\nCompiling FAST without openmp\n');
    make('openmp',false,'fsparseonly',1);

end

cd ../

end
