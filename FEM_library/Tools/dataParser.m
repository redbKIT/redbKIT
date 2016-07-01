function [ DATA ] = dataParser( DATA )
%DATAPARSER input parser

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch> 


DATA = parserLinearSolverOptions( DATA );
DATA = parserPreconditionerOptions( DATA );
DATA = parserNonLinearSolverOptions( DATA );

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ DATA ] = parserLinearSolverOptions( DATA )

% Set Default options if not provided

if ~isfield(DATA, 'LinearSolver')
    DATA.LinearSolver = [];
end

if ~isfield(DATA.LinearSolver,'type')
    DATA.LinearSolver.type                =  'backslash';
end

if ~isfield(DATA.LinearSolver,'gmres_verbosity')
    DATA.LinearSolver.gmres_verbosity     =  1;
end

if ~isfield(DATA.LinearSolver,'tol')
    DATA.LinearSolver.tol                 =  1e-8;
end

if ~isfield(DATA.LinearSolver,'maxit')
    DATA.LinearSolver.maxit               =  150;
end

if ~isfield(DATA.LinearSolver,'mumps_reordering')
    DATA.LinearSolver.mumps_reordering    =  0;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ DATA ] = parserPreconditionerOptions( DATA )

if ~isfield(DATA, 'Preconditioner')
    DATA.Preconditioner = [];
end

if ~isfield(DATA.Preconditioner,'type')
    DATA.Preconditioner.type                        =  'None';
end

if ~isfield(DATA.Preconditioner,'local_solver')
    DATA.Preconditioner.local_solver                =  'matlab_lu';
end

if ~isfield(DATA.Preconditioner,'overlap_level')
    DATA.Preconditioner.overlap_level               =  1;
end

if ~isfield(DATA.Preconditioner,'mumps_reordering')
    DATA.Preconditioner.mumps_reordering            =  0;
end

if ~isfield(DATA.Preconditioner,'num_subdomains')
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        poolsize = 0;
    else
        poolsize = poolobj.NumWorkers;
    end
    DATA.Preconditioner.num_subdomains              =  poolsize;
end

if ~isfield(DATA.Preconditioner,'ILU_type')
    DATA.Preconditioner.ILU_type            =  'ilutp';
end

if ~isfield(DATA.Preconditioner,'ILU_droptol')
    DATA.Preconditioner.ILU_droptol         =  1e-2;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ DATA ] = parserNonLinearSolverOptions( DATA )

if ~isfield(DATA, 'NonLinearSolver')
    DATA.NonLinearSolver = [];
end

if ~isfield(DATA.NonLinearSolver,'type')
    DATA.NonLinearSolver.type          =  'Newton';
end

if ~isfield(DATA.NonLinearSolver,'tol')
    DATA.NonLinearSolver.tol           =  1e-5;
end

if ~isfield(DATA.NonLinearSolver,'maxit')
    DATA.NonLinearSolver.maxit         =  15;
end

if ~isfield(DATA.NonLinearSolver,'backtrackIter')
    DATA.NonLinearSolver.backtrackIter         =  0;
end

if ~isfield(DATA.NonLinearSolver,'backtrackFactor')
    DATA.NonLinearSolver.backtrackFactor       =  1.0;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
