%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.SaturateOnIntegerOverflow = false;
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = false;
cfg.ExtrinsicCalls = false;

%% Define argument types for entry-point 'STR_assembler2D_quad'.
ARGS = cell(1,1);
ARGS{1} = cell(11,1);
ARGS{1}{1} = coder.typeof(0,[1 1],[0 0]);%dim
ARGS{1}{2} = coder.typeof(0,[1 1],[0 0]);%Material_Model
ARGS{1}{3} = coder.typeof(0,[1  50],[0 1]);%material_param
ARGS{1}{4} = coder.typeof(0,[Inf 1],[1 0]);%Uh
ARGS{1}{5} = coder.typeof(0,[40 inf],[1 1]);%elements
ARGS{1}{6} = coder.typeof(0,[1 1],[0 0]);%nln
ARGS{1}{7} = coder.typeof(0,[1 100],[0 1]);%w
ARGS{1}{8} = coder.typeof(0,[inf 3 3],[1 1 1]);%invjac
ARGS{1}{9} = coder.typeof(0,[1 inf],[0 1]);%detjac
ARGS{1}{10} = coder.typeof(0,[40 100],[1 1]);%phi
ARGS{1}{11} = coder.typeof(0,[40 100 3],[1 1 1]);%dphi

%% Invoke MATLAB Coder.
cfg.PostCodeGenCommand = sprintf('buildInfo.addLinkFlags(''-fopenmp'')');
codegen -config cfg CSM_assembler_M_omp -args ARGS{1}

