function cmd = bfsl_ica_aroma(infile, TR, varargin)
    % BFSL_ICA_AROMA Run ICA-AROMA using FSL
    %
    %  USAGE: bfsl_ica_aroma(infile, TR, outdir)
    % __________________________________________________________________________
    %  INPUTS
    %   infile:     4D timeseries to denoise
    %   TR:         in seconds [optional]
    %   outdir:     output directory [default = same as in]
    %   dentype:    Type of denoising strategy (default is nonaggr):
    %               no: only classification, no denoising
    %               nonaggr: non-aggressive (partial component regression)
    %               aggr: aggressive (full component regression)
    %               both: aggressive and non-aggressive (two outputs)
    %   affmat      File name of the mat-file describing the affine registration
    %               (e.g. FSL FLIRT) of the functional data to structural space
    %               (.mat file)
    %   warpfile    name of the warp-file describing the non-linear registration
    %               (e.g. FSL FNIRT) of the structural data to MNI152 space (.nii.gz)
    %
    % USAGE: ICA_AROMA.py [-h] -o OUTDIR [-i INFILE] [-mc MC] [-a AFFMAT] [-w WARP]
    %                     [-m MASK] [-f INFEAT] [-tr TR] [-den DENTYPE] [-md MELDIR]
    %                     [-dim DIM]
    %
    % EXAMPLE: ICA_AROMA.py -in filtered_func_data.nii.gz -out ICA_AROMA -mc rest_mcf.par
    %                       -m mask_aroma.nii.gz -affmat func2highres.mat -warp
    %                       highres2standard_warp.nii.gz -md filtered_func_data.ica

    % ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
    %   Created:  2015-03-24
    %   Email:    spunt@caltech.edu
    % __________________________________________________________________________
    def = { ...
        'outdir',        [],                  ...
        'dentype',       'nonaggr',                 ...
        'affmat',        [],                 ...
        'warpfile',      [],                  ...
         };
    vals = setargs(def, varargin);
    if nargin < 2, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
    if iscell(infile), infile = char(infile); end
    [p,n,e] = fileparts(infile);
    if ~strcmp(e, '.gz'), pigz(infile); infile = strcat(infile, '.gz'); end
    if isempty(outdir), outdir = fileparts(infile); end
    rpfile      = char(files(fullfile(fileparts(infile), 'rp*txt')));
    if isempty(rpfile), disp('Motion Correction file not found!'); return; end
    % | - Configure Path
    gitdir      = fullfile(getenv('HOME'), 'Github');
    aromadir    = fullfile(gitdir, 'thirdparty-matlab', 'ICA-AROMA');
    icaaroma    = fullfile(aromadir, 'ICA_AROMA.py');
    if ~exist(icaaroma, 'file'), fprintf('\n\nPATH IS INVALID: %s\n', icaaroma); return; end
    % | - Construct Command
    cmd = sprintf('python %s -i %s -o %s -mc %s -tr %d -den %s', icaaroma, infile, outdir, rpfile, TR, dentype);
    if ~isempty(affmat)
        if iscell(affmat), affmat = char(affmat); end
        cmd = sprintf('%s -a %s', cmd, affmat);
    end
    if ~isempty(warpfile)
        if iscell(warpfile), warpfile = char(warpfile); end
        [p,n,e] = fileparts(warpfile);
        if ~strcmp(e, '.gz'), pigz(warpfile); warpfile = strcat(warpfile, '.gz'); end
        cmd = sprintf('%s -w %s', cmd, warpfile);
    end
    if nargout==0, system(cmd); end
function argstruct = setargs(defaultargs, varargs)
    if nargin < 1, mfile_showhelp; return; end
    if nargin < 2, varargs = []; end
    defaultargs = reshape(defaultargs, 2, length(defaultargs)/2)';
    if ~isempty(varargs)
        if mod(length(varargs), 2)
            error('Optional inputs must be entered as "''Name'', Value" pairs, e.g., myfunction(''arg1'', val1, ''arg2'', val2)');
        end
        arg = reshape(varargs, 2, length(varargs)/2)';
        for i = 1:size(arg,1)
           idx = strncmpi(defaultargs(:,1), arg{i,1}, length(arg{i,1}));
           if sum(idx) > 1
               error(['Input "%s" matches multiple valid inputs:' repmat('  %s', 1, sum(idx))], arg{i,1}, defaultargs{idx, 1});
           elseif ~any(idx)
               error('Input "%s" does not match a valid input.', arg{i,1});
           else
               defaultargs{idx,2} = arg{i,2};
           end
        end
    end
    for i = 1:size(defaultargs,1), assignin('caller', defaultargs{i,1}, defaultargs{i,2}); end
    if nargout>0, argstruct = cell2struct(defaultargs(:,2), defaultargs(:,1)); end
function mfile_showhelp(varargin)
    % MFILE_SHOWHELP
    ST = dbstack('-completenames');
    if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
    eval(sprintf('help %s', ST(2).file));

