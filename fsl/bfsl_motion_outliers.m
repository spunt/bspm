function cmd = bfsl_motion_outliers(epi, TR, varargin)
    %  USAGE: bfsl_motion_outliers(epi, varargin);
    %
    %  metric: 'refrms', 'dvars', 'refmse', 'fd', 'fdrms'
    % __________________________________________________________________________
    % Usage: fsl_motion_outliers -i <input 4D image> -o <output confound file> [options]
    %
    % Options:
    %      -m <mask image>      use supplied mask image for calculating metric
    %      -s <filename>        save metric values (e.g. DVARS) as text into specified file
    %      -p <filename>        save metric values (e.g. DVARS) as a graphical plot (png format)
    %      -t <path>            [Optional] Path to the location where temporary files should be created. Defaults to /tmp
    %      --refrms             use RMS intensity difference to reference volume as metric [default metric]
    %      --dvars              use DVARS as metric
    %      --refmse             Mean Square Error version of --refrms
    %      --fd                 use FD (framewise displacement) as metric
    %      --fdrms              use FD with RMS matrix calculation as metric
    %      --thresh=<val>       specify absolute threshold value (otherwise use box-plot cutoff = P75 + 1.5*IQR)
    %      --nomoco             do not run motion correction (assumed already done)
    %      --dummy=<val>        number of dummy scans to delete (before running anything and creating EVs)
    %      -v                   verbose mode
    %
    def = {
    'include_rp', 1,             ...
    'makeplot',   0,             ...
    'maskfile',   [],            ...
    'metric',     'refrms',      ...
    'nomoco',     true,          ...
    'prefix',     'fslbadscan',  ...
    'saveraw',    0,             ...
    'thresh',     'auto',        ...
    };
    vals = setargs(def, varargin);
    if nargin < 1, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
    if iscell(epi), epi = char(epi); end
    [epidir,n,e] = fileparts(epi);
    if ~strcmp(e, '.gz'), pigz(epi); epi = strcat(epi, '.gz'); end

    % | create tmp output filename
    tmpfile = fullfile(epidir, 'tmp.txt');

    % | - Construct Command
    cmd1 = sprintf('fsl_motion_outliers -i %s -o %s', epi, tmpfile);
    if ~isempty(maskfile)
        if iscell(maskfile), maskfile = char(maskfile); end
        [p,n,e] = fileparts(maskfile);
        if ~strcmp(e, '.gz'), pigz(maskfile); maskfile = strcat(maskfile, '.gz'); end
        cmd1 = sprintf('%s -m %s', cmd1, maskfile);
    end
    cmd2 = sprintf('--%s', metric);
    if nomoco, cmd2 = sprintf('%s --nomoco', cmd2); end
    if ~strcmpi(thresh, 'auto'), cmd2 = sprintf('%s --thresh=%2.2f', cmd2, thresh); end
    cmd = sprintf('%s %s', cmd1, cmd2);
    system(cmd);
    if ~exist(tmpfile), fprintf('\n\nDAMN! SOMETHING WENT WRONG...\n\n'); return; end;
    scrubmat = load(tmpfile);
    [nvol, nbad] = size(scrubmat);
    pctbad = round(100*(nbad/nvol));
    fprintf('\n | - %d (%d%%) scans identified as bad\n', nbad, pctbad);

    % | motion parameters
    if include_rp
        rpfile      = char(files(fullfile(fileparts(epi), 'rp*txt')));
        if isempty(rpfile), disp('Motion Correction file not found!'); return; end
        rp                  = load(rpfile);
        rp(:,4:6)           = rp(:,4:6)*57.2958;
        scrubmat = [rp scrubmat];
    end
    outfile     = sprintf('%s_%s_%dpctbad_%s.txt', prefix, metric, pctbad, strtrim(datestr(now,'mmm_DD_YYYY')));

    % | save as a new text file
    save(fullfile(epidir, outfile), 'scrubmat', '-ascii');
    fprintf('\n | - Nuisance regressors written to:  %s\n\n', outfile);
    delete(tmpfile);

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

