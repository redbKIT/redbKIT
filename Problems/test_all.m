function test_all
%TEST_ALL launch all the tests and create a log file

%   This file is part of redbKIT.
%   Copyright (c) 2015, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 


fid = fopen('test_log.txt','w');

TestFolder = pwd;

%% Write file header
t = now;
c = datevec ( t );
s = datestr ( c, 0 );
fprintf(fid, 'Test logfile: %s\n', s );
s = version;
fprintf(fid, '\nMatlab Version: %s\n', s );
s = computer;
fprintf(fid, 'Platform: %s\n\n', s );

%% Start Testing

% Test
cd FEM_Test_2DLaplacian/
try
    testname = 'FEM_Test_2DLaplacian_P1';
    test('P1','Dirichlet_data');    
    test('P1','DirichletNeumann_data');
    test('P1','DirichletRobin_data');
    close all;
    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_Test_2DLaplacian_P2';
    test('P2','Dirichlet_data');    
    test('P2','DirichletNeumann_data');
    test('P2','DirichletRobin_data');
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd FEM_Test_3DLaplacian/
try
    testname = 'FEM_Test_3DLaplacian_P1';
    test('P1','Dirichlet_data');    
    close all;
    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_Test_3DLaplacian_P2';
    test('P2','Dirichlet_data');    
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd FEM_Test_ADRt_2D/
try
    testname = 'FEM_Test_ADRt_2D_P1';
    Solve('P1');    
    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_Test_ADRt_2D_P2';
    Solve('P2');    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd FEM_Test_ADRt_3D/
try
    testname = 'FEM_Test_ADRt_3D_P1';
    Solve('P1');    
    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd FEM_CSM_Test2D/
try
    testname = 'FEM_CSM_Test2D_P1';
    test('P1');    
    close all;
    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CSM_Test2D_P2';
    test('P2');    
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd FEM_CSM_Test3D/
try
    testname = 'FEM_CSM_Test3D_P1';
    test('P1');    
    close all;
    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CSM_Test3D_P2';
    test('P2');    
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)


% Test
cd FEM_CSMt_Test2D/
try
    testname = 'FEM_CSMt_Test2D_LinearP1';
    Solve('Linear', 3, 'P1');    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CSMt_Test2D_SVKP1';
    Solve('StVenantKirchhoff', 3, 'P1');    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CSMt_Test2D_LinearP2';
    Solve('Linear', 3, 'P2');    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CSMt_Test2D_SVKP2';
    Solve('StVenantKirchhoff', 3, 'P2');    
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test CFD 2D
cd FEM_CFD_Dfg2D/
try
    testname = 'FEM_CFD_Dfg2D_P2ImplicitNoStab';
    test('P2', 'implicit', 'None');
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CFD_Dfg2D_P2SemiImplicitNoStab';
    test('P2', 'semi-implicit', 'None');
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CFD_Dfg2D_P1SemiImplicitSUPG';
    test('P1', 'semi-implicit', 'SUPG');
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CFD_Dfg2D_P1ImplicitSUPG';
    test('P1', 'implicit', 'SUPG');
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test CFD 3D
cd FEM_CFD_Dfg3D/
try
    testname = 'FEM_CFD_Dfg3D_P2ImplicitNoStab';
    test('B1', 'implicit', 'None');
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CFD_Dfg3D_P2SemiImplicitNoStab';
    test('P2', 'semi-implicit', 'None');
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CFD_Dfg3D_P1SemiImplicitSUPG';
    test('P1', 'semi-implicit', 'SUPG');
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
try
    testname = 'FEM_CFD_Dfg3D_P1ImplicitSUPG';
    test('P1', 'implicit', 'SUPG');
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd FEM_TestMetis/
try
    testname = 'FEM_TestMetis';
    test(2);
    test(6);
    test(12);
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd Test_EIM_DEIM/
try
    testname = 'Test_EIM_DEIM';
    test;
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd RB_AffineDevice/
try
    testname = 'RB_AffineDevice';
    test;
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd RB_EIM_Gaussian/
try
    testname = 'RB_EIM_Gaussian';
    test;
    close all;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

% Test
cd RB_CSM_ShearCube/
try
    testname = 'RB_CSM_ShearCube';
    test;
    print_test_passed(fid,testname);
    
catch err
    print_error_toFile(fid,err,testname);
end
cd(TestFolder)

%% End Testing

t = now;
c = datevec ( t );
s = datestr ( c, 0 );
fprintf(fid,'----------------------------------------------------\n');
fprintf(fid, '\nFinished: %s\n', s );

fclose(fid);

fprintf('\n----------------------------------------------------\n');
fprintf('----------------------------------------------------\n');
type('test_log.txt')

end



function print_error_toFile(fid,err,testname)

fprintf(fid,'-------------------------------\n');
fprintf(fid,['Error in ',testname,'\n\n']);

fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
fprintf(fid,'\n');

end

function print_test_passed(fid,testname)
fprintf(fid,'----------------------------------------------------\n');
fprintf(fid,'%38s | passed |\n',testname);

end
