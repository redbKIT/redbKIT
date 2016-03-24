%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script generated from project 'ADR_assembly.prj' on 10-Dec-2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create configuration object of class 'coder.MexCodeConfig'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.SaturateOnIntegerOverflow = false;
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = false;
cfg.ExtrinsicCalls = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define argument types for entry-point 'ADR_mex_assembler'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ARGS = cell(1,1);
ARGS{1} = cell(26,1);
ARGS{1}{1} = coder.typeof('X',[1 50],[0 1]);%OPERATOR
ARGS{1}{2} = coder.typeof(0,[1 2],[0 1]);%TC_d
ARGS{1}{3} = coder.typeof(0);%TC_t
ARGS{1}{4} = coder.typeof(0,[20 inf],[1 1]);%elements
ARGS{1}{5} = coder.typeof(0);%nln
ARGS{1}{6} = coder.typeof(0,[inf    120],[1 1]);%mu
ARGS{1}{7} = coder.typeof(0,[inf    120],[1 1]);%bx
ARGS{1}{8} = coder.typeof(0,[inf    120],[1 1]);%by
ARGS{1}{9} = coder.typeof(0,[inf    120],[1 1]);%bz
ARGS{1}{10} = coder.typeof(0,[inf    120],[1 1]);%si
ARGS{1}{11} = coder.typeof(0,[inf    120],[1 1]);%f
ARGS{1}{12} = coder.typeof(0,[1 120],[0 1]);%w
ARGS{1}{13} = coder.typeof(0,[1 inf],[0 1]);%dcdx
ARGS{1}{14} = coder.typeof(0,[1 inf],[0 1]);%dcdy
ARGS{1}{15} = coder.typeof(0,[1 inf],[0 1]);%dcdz
ARGS{1}{16} = coder.typeof(0,[1 inf],[0 1]);%dedx
ARGS{1}{17} = coder.typeof(0,[1 inf],[0 1]);%dedy
ARGS{1}{18} = coder.typeof(0,[1 inf],[0 1]);%dedz
ARGS{1}{19} = coder.typeof(0,[1 inf],[0 1]);%dtdx
ARGS{1}{20} = coder.typeof(0,[1 inf],[0 1]);%dtdy
ARGS{1}{21} = coder.typeof(0,[1 inf],[0 1]);%dtdz
ARGS{1}{22} = coder.typeof(0,[20 120],[1 1]);%phi
ARGS{1}{23} = coder.typeof(0,[20 120],[1 1]);%dcsiphi
ARGS{1}{24} = coder.typeof(0,[20 120],[1 1]);%detaphi
ARGS{1}{25} = coder.typeof(0,[20 120],[1 1]);%dtauphi
ARGS{1}{26} = coder.typeof(0,[1 inf],[0 1]);%detjac

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Invoke MATLAB Coder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.PostCodeGenCommand = sprintf('buildInfo.addLinkFlags(''-fopenmp'')');
codegen -config cfg ADR_mex_assembler3D.m -args ARGS{1}
