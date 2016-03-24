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
ARGS{1} = cell(19,1);
ARGS{1}{1} = coder.typeof('X',[1 50],[0 1]);
ARGS{1}{2} = coder.typeof(0,[1 2],[0 1]);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0,[20 inf],[1 1]);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0,[inf    60],[1 1]);
ARGS{1}{7} = coder.typeof(0,[inf    60],[1 1]);
ARGS{1}{8} = coder.typeof(0,[inf    60],[1 1]);
ARGS{1}{9} = coder.typeof(0,[inf    60],[1 1]);
ARGS{1}{10} = coder.typeof(0,[inf    60],[1 1]);
ARGS{1}{11} = coder.typeof(0,[1 60],[0 1]);
ARGS{1}{12} = coder.typeof(0,[1 inf],[0 1]);
ARGS{1}{13} = coder.typeof(0,[1 inf],[0 1]);
ARGS{1}{14} = coder.typeof(0,[1 inf],[0 1]);
ARGS{1}{15} = coder.typeof(0,[1 inf],[0 1]);
ARGS{1}{16} = coder.typeof(0,[20 60],[1 1]);
ARGS{1}{17} = coder.typeof(0,[20 60],[1 1]);
ARGS{1}{18} = coder.typeof(0,[20 60],[1 1]);
ARGS{1}{19} = coder.typeof(0,[1 inf],[0 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Invoke MATLAB Coder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.PostCodeGenCommand = sprintf('buildInfo.addLinkFlags(''-fopenmp'')');
codegen -config cfg ADR_mex_assembler -args ARGS{1}
