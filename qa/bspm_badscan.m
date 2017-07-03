function tsinfo = bspm_badscan(epi, varargin)
    % BSPM_BADSCAN Identifies bad scans and saves a nuisance regressor file
    %
    %   USAGE: tsinfo = bspm_badscan(epi, varargin)
    %
    %  ARGUMENTS
    %    epi = cell array of EPI filenames
    %   'dvars_thresh',     2.5,
    %   'framewise_thresh', 0.5,
    %   'include_rp',       1,
    %   'prefix',           'badscan',
    %   'maskfile',         []
    %
    %
    % Bob Spunt, Caltech
    % Based on a script by Donald McLaren, Ph.D.
    % 2012_05_01 -- Created
    % 2012_06_04 -- Added to FUNC + Presented in SCAN Lab Workshop
    % 2013_03_08 -- Added option to MASK data prior to computing bad scans
    % ========================================================================%
    DOONEOUTZSCORE = false;
    def = { 'dvars_thresh',     2.5,  ...
            'framewise_thresh', 0.5,  ...
            'include_rp',       1,    ...
            'makeplot',         0,    ...
            'prefix',           'badscan', ...
            'maskfile',         []};
    vals = setargs(def, varargin);
    if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end

    % | PATH TO BRAMILA TOOLS
    bramiladir = fullfile(getenv('HOME'), 'Github', 'bspm', 'thirdparty', 'bramila');
    addpath(bramiladir)

    % | get data and apply implicit masking
    if ischar(epi), epi = cellstr(epi); end
    fprintf('\n | - Loading Data');
    maskthresh  = 0.8;
    [v,h]       = bnii_read(epi, 'reshape', 1);
    [nvox nvol] = size(v);
    if isempty(maskfile)
        v(v < repmat(mean(v)*maskthresh, nvox, 1)) = NaN;
    else
        m = bspm_reslice(maskfile, strcat(epi, ',1'), 0, 1);
        v(m(:)==0,:) = NaN;
    end
    cfg.vol     = reshape(v, h.dim(2:5));

    % | create output filename
    outfile     = sprintf('%s_dvars%dframewise%d_%s.txt', prefix, dvars_thresh*100, framewise_thresh*100, strtrim(datestr(now,'mmm_DD_YYYY')));
    epidir      = fileparts(epi{1});

    % | motion parameters
    rp                  = load(char(files(fullfile(epidir, 'rp_*txt'))));
    rp(:,4:6)           = rp(:,4:6)*57.2958;
    framewise           = zeros(nvol,1);
    framewise(2:end)    = max(abs(diff(rp)), [], 2);

    % | use BRAMILA tools to compute dvars and framewise displacement
    fprintf('\n | - Identifying Scans with DVARS > %2.2f and framewise displacement > %2.2f\n', dvars_thresh, framewise_thresh);
    cfg.plot            = makeplot;
    dvars               = bramila_dvars(cfg);

    % | get indices of bad timepoints
    if DOONEOUTZSCORE
        zdvars          = oneoutzscore(dvars(2:end), 1);
    else
        zdvars          = zscore(dvars(2:end), 1);
    end
    badidx          = zeros(nvol, 2);
    badidx(2:end,1) = zdvars > dvars_thresh;
    badidx(:,2)     = framewise > framewise_thresh;
    badidx          = find(any(badidx, 2));
    nbad            = length(badidx);
    tsinfo.pctbad   = round(100*(nbad/nvol));
    fprintf('\n | - %d (%d%%) scans identified as bad', nbad, tsinfo.pctbad);

    % | create bad volume regressor
    scrubmat = zeros(nvol, nbad);
    for i = 1:length(badidx), scrubmat(badidx(i),i) = 1; end
    if include_rp==1, scrubmat = [rp scrubmat]; end

    % | save as a new text file
    save(fullfile(epidir, outfile), 'scrubmat', '-ascii');
    fprintf('\n | - Nuisance regressors written to:  %s\n\n', outfile);

    if nargout>0
        tsinfo.epidir       = epidir;
        tsinfo.maxmotion    = max(rp) - min(rp);
        tsinfo.dvars        = dvars;
        tsinfo.framewise    = framewise;
        tsinfo.scrubmat     = scrubmat;
    end
    rmpath(bramiladir)
function zx = oneoutzscore(x, returnas)
    % ONEOUTZSCORE Perform columnwise leave-one-out zscoring
    %
    % USAGE: zx = oneoutzscore(x, returnas)
    %
    %   returnas: 0, signed values (default); 1, absolute values
    %
    if nargin<1, disp('USAGE: zx = oneoutzscore(x, returnas)'); return; end
    if nargin<2, returnas = 0; end
    if size(x,1)==1, x=x'; end
    zx              = x;
    [nrow, ncol]    = size(x);
    for c = 1:ncol
        cin         = repmat(x(:,c), 1, nrow);
        theoneout   = cin(logical(eye(nrow)))';
        theleftin   = reshape(cin(logical(~eye(nrow))),nrow-1,nrow);
        cz          = (theoneout-nanmean(theleftin))./nanstd(theleftin);
        zx(:,c)     = cz';
    end
    if returnas, zx = abs(zx); end




