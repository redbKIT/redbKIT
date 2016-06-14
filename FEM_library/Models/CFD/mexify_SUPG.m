%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

%==========================================================================
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
ARGS{1} = cell(14,1);
ARGS{1}{1} = coder.typeof(0,[1 1],[0 0]);%dim
ARGS{1}{2} = coder.typeof(0,[40 Inf],[1 1]);%elements
ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]);%detjac
ARGS{1}{4} = coder.typeof(0,[Inf 3 3],[1 1 1]);%invjac
ARGS{1}{5} = coder.typeof(0,[1 100],[0 1]);%w
ARGS{1}{6} = coder.typeof(0,[40 100],[1 1]);%phi
ARGS{1}{7} = coder.typeof(0,[40 100 3],[1 1 1]);%dphiref
ARGS{1}{8} = coder.typeof(0,[1 1],[0 0]);%nln
ARGS{1}{9} = coder.typeof(0,[Inf 1],[1 0]);%uh
ARGS{1}{10} = coder.typeof(0,[Inf 1],[1 0]);%un
ARGS{1}{11} = coder.typeof(0,[1 1],[0 0]);%density
ARGS{1}{12} = coder.typeof(0,[1 1],[0 0]);%viscosity
ARGS{1}{13} = coder.typeof(0,[1 1],[0 0]);%dt
ARGS{1}{14} = coder.typeof(0,[1 1],[0 0]);%alpha

%% Invoke MATLAB Coder.
cfg.PostCodeGenCommand = sprintf('buildInfo.addLinkFlags(''-fopenmp'')');
codegen -config cfg CFD_Assemble_SUPG_SemiImpl -args ARGS{1}

%==========================================================================
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
ARGS{1} = cell(14,1);
ARGS{1}{1} = coder.typeof(0,[1 1],[0 0]);%dim
ARGS{1}{2} = coder.typeof(0,[40 Inf],[1 1]);%elements
ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]);%detjac
ARGS{1}{4} = coder.typeof(0,[Inf 3 3],[1 1 1]);%invjac
ARGS{1}{5} = coder.typeof(0,[1 100],[0 1]);%w
ARGS{1}{6} = coder.typeof(0,[40 100],[1 1]);%phi
ARGS{1}{7} = coder.typeof(0,[40 100 3],[1 1 1]);%dphiref
ARGS{1}{8} = coder.typeof(0,[1 1],[0 0]);%nln
ARGS{1}{9} = coder.typeof(0,[Inf 1],[1 0]);%U_k
ARGS{1}{10} = coder.typeof(0,[Inf 1],[1 0]);%v_n
ARGS{1}{11} = coder.typeof(0,[1 1],[0 0]);%density
ARGS{1}{12} = coder.typeof(0,[1 1],[0 0]);%viscosity
ARGS{1}{13} = coder.typeof(0,[1 1],[0 0]);%dt
ARGS{1}{14} = coder.typeof(0,[1 1],[0 0]);%alpha

%% Invoke MATLAB Coder.
cfg.PostCodeGenCommand = sprintf('buildInfo.addLinkFlags(''-fopenmp'')');
codegen -config cfg CFD_Assemble_SUPG_Impl -args ARGS{1}

