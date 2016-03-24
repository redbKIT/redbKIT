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
ARGS{1} = cell(3,1);
ARGS{1}{1} = coder.typeof(0);%dim
ARGS{1}{2} = coder.typeof(0,[3  inf],[1 1]);%vertices
ARGS{1}{3} = coder.typeof(0,[20 inf],[1 1]);%elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Invoke MATLAB Coder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.PostCodeGenCommand = sprintf('buildInfo.addLinkFlags(''-fopenmp'')');
codegen -config cfg geotrasf.m -args ARGS{1}
