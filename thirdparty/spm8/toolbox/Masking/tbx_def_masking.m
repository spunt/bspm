function varargout = tbx_def_masking(defstr, defval)
% MATLABBATCH Defaults file for toolbox 'Mask Creation'
% This code has been automatically generated.

persistent defaults;
if isempty(defaults)
    defaults.makeavg.avgexpr = 'mean(X)';
    defaults.makeavg.outname = 'average.nii';
    defaults.makeavg.outdir = '';
    defaults.optthr.optfunc = '@opt_thr_corr';
    defaults.optthr.outname = 'average_optthr.nii';
    defaults.optthr.outdir = '';
end
if nargin == 1
    [un varargout{1}] = evalc(['defaults.' defstr ';']);
else
    evalc(['defaults.' defstr ' = defval;']);
end
