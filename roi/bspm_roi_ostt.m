function varargout = bspm_roi_ostt(level1dirs, roifiles, conidx, varargin)
% BSPM_ROI_OSTT Test level 1 contrasts
%
%  USAGE: bspm_roi_ostt(level1dirs, roifiles, conidx, varargin)
% ________________________________________________________________________________________
%  INPUTS
%	level1dirs:  
%	roifiles:  
%	conidx:  
%

% ----------------------------- Copyright (C) 2016 Bob Spunt -----------------------------
%	Created:  2016-05-16
%	Email:     spunt@caltech.edu
% ________________________________________________________________________________________
def = { ...
'contrasts',            []                      , ...
'rmsuboutlier',         0                       , ...
'roinames',             ''                      , ...
'dofdr',                1                       , ...
'condirpat',            ''                      , ...
'tag',                  ''                      , ...
'savetable',            1                       , ...
'conpat',               'con*nii*'              , ...
'templatefile'          ''  ...
	};
vals = setargs(def, varargin);
if nargin < 3, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end

% | CHECK REQUIRED INPUTS
% | ======================================================================================
if ischar(level1dirs), level1dirs = cellstr(level1dirs); end
if ischar(roifiles), roifiles = cellstr(roifiles); end
nlevel1 = length(level1dirs); 
nroi    = length(roifiles);
if any(~cellfun(@exist, level1dirs, repmat({'dir'}, nlevel1, 1)))
    printmsg('One or more level 1 directory does not exist!'); 
    return; 
end
if any(~cellfun(@exist, roifiles, repmat({'file'}, nroi, 1)))
    printmsg('One or more ROI files does not exist!'); 
    return; 
end

% | CHECK VARARGIN
% | ======================================================================================
[nrow, ncol] = size(conidx);
ncond = numel(conidx);
if isempty(roinames), [~,roinames] = cellfileparts(roifiles); end
if length(roinames)~=nroi, printmsg('# of roinames must match # of rois'); return; end
if isempty(contrasts), contrasts = eye(ncond); end
if size(contrasts, 2) > ncond, printmsg('too many columns in contrasts variable!'); return; end
ncontrast       = size(contrasts,1);
% | GET CONTRAST FILES
% | ======================================================================================
conpat = fullfile(level1dirs, conpat);
allcon = cellfun(@files, conpat, 'Unif', false); 
allcon = cellfun(@(x) x(conidx(:)), allcon, 'Unif', false);
allcon = vertcat(allcon{:});
if ncond < max(conidx(:))
    tmpidx              = zeros(ncond, 2);
    tmpidx(:,1)         = 1:ncond;
    tmpidx(:,2)         = conidx(:);
    tmpidx              = sortrows(tmpidx, 2);
    conidx(tmpidx(:,1)) = 1:ncond;
end

% | READ IMAGE VOLUME FILES
% | ======================================================================================
[CON, CONHDR] = bspm_read_vol(allcon, 'reshape');
[ROI, ROIHDR] = bspm_read_vol(roifiles, 'reshape');
condname = regexprep({CONHDR(1:ncond).descrip}', '.+\d:\s', '');
contrastnames = bspm_conweights2names(contrasts, condname);
contrastnames = regexprep(contrastnames, '_-_', ' > ');
roinames = regexprep(roinames, '\W', ' '); 
brainname = [];
for c = 1:ncontrast, brainname = [brainname; strcat(upper(roinames), {' | '}, contrastnames(c))]; end

% | LOOP OVER ROIS
% | ======================================================================================
allrow = [];
colsep = repmat({''}, nroi, 1);
for i = 1:nroi
    cdata       = reshape(nanmean(CON(ROI(:,i)>0,:)), ncond, nlevel1)';
    roirow      = [];
    roirow{1}   = roinames{i};
    for c = 1:size(contrasts,1)
        contrast     = contrasts(c,:);
        cmat         = repmat(contrast, size(cdata,1), 1);
        contrastdata = cdata(:,abs(contrast)>0);
        y = sum(contrastdata.*cmat(:,abs(contrast)>0), 2);
        if rmsuboutlier, y = bob_outlier2nan(y, rmsuboutlier); end      
        out(c)       = bob_ttest(y(~isnan(y)), 1, contrastnames{c});
        roirow       = [roirow {''} num2cell([out(c).tstat out(c).P out(c).CI'])];
    end
    allrow = [allrow; roirow];
end
if dofdr
    pidx    = 4:5:(size(allrow,2)); 
    pidx    = pidx(1:size(contrasts,1))
    pval    = allrow(:,pidx);
    pval    = cell2mat(pval);
    fdrpval = pval;
    for i = 1:size(contrasts,1)
        [h, crit_p(i), fdrpval(:,i)] = fdr_bh(pval(:,i), .05, 'dep', 'yes');
    end
    allrow(:,pidx) = num2cell(fdrpval);
end

% | SAVE
% | ======================================================================================
if savetable
    if isempty(templatefile)
        templatefile = fullfile(fileparts(mfilename('fullpath')), 'TEMPLATE_ROI_OSTT.xlsx');
    else
        if iscell(templatefile), templatefile = char(templatefile); end
        if ~exist(templatefile, 'file')
            printmsg('Custom template not found on path! Using default...'); 
            templatefile = fullfile(fileparts(mfilename('fullpath')), 'TEMPLATE_ROI_OSTT.xlsx');
        end
    end
    hdr                   = cell(4, size(allrow, 2)); 
    hdr(2,3:5:end)        = contrastnames;
    hdr{3, 1}             = 'ROI'; 
    if dofdr
        vhdr = {'t-stat' 'pFDR' '95% CI' '' ''};
    else
        vhdr = {'t-stat' 'pFDR' '95% CI' '' ''};
    end
    vhdr                  = repmat(vhdr, 1, ncontrast); 
    hdr(3,3:end)          = vhdr(1:end-1);
    TABLE                 = [hdr; allrow]; 
    outstr                = {sprintf('Table_ROI_Analysis_%s', tag) sprintf('_SubOut-%dSD', rmsuboutlier)};
    outname               = strcat(outstr{find([1 rmsuboutlier])});
    outname               = regexprep(outname, '_$', '');
    outname               = fullfile(pwd, strcat(outname, '.xlsx'));
    
    copyfile(templatefile, outname);
    xlwrite(outname, TABLE);
end

end
% ========================================================================================
%
% ------------------------------------- SUBFUNCTIONS -------------------------------------
%
% ========================================================================================
function mfile_showhelp(varargin)
% MFILE_SHOWHELP
ST = dbstack('-completenames');
if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
eval(sprintf('help %s', ST(2).file));
end
function argstruct = setargs(defaults, optargs)
% SETARGS Name/value parsing and assignment of varargin with default values
if nargin < 1, mfile_showhelp; return; end
if nargin < 2, optargs = []; end
defaults = reshape(defaults, 2, length(defaults)/2)'; 
if ~isempty(optargs)
    if mod(length(optargs), 2)
        error('Optional inputs must be entered as Name, Value pairs, e.g., myfunction(''name'', value)'); 
    end
    arg = reshape(optargs, 2, length(optargs)/2)';
    for i = 1:size(arg,1)
       idx = strncmpi(defaults(:,1), arg{i,1}, length(arg{i,1}));
       if sum(idx) > 1
           error(['Input "%s" matches multiple valid inputs:' repmat('  %s', 1, sum(idx))], arg{i,1}, defaults{idx, 1});
       elseif ~any(idx)
           error('Input "%s" does not match a valid input.', arg{i,1});
       else
           defaults{idx,2} = arg{i,2};
       end  
    end
end
for i = 1:size(defaults,1), assignin('caller', defaults{i,1}, defaults{i,2}); end
if nargout>0, argstruct = cell2struct(defaults(:,2), defaults(:,1)); end
end
