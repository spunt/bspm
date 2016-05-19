function varargout = bspm_plot_contrast(level1dirs, roifiles, conidx, varargin)
% BSPM_PLOT_CONTRAST Plot level 1 contrasts
%
%  USAGE: bspm_plot_contrast(level1dirs, roifiles, conidx, varargin)
% ________________________________________________________________________________________
%  INPUTS
%	level1dirs:  
%	roifiles:  
%	conidx:  
%
% ________________________________________________________________________________________
%  VARARGIN
%	vararg1:  
%	vararg2:  
%	vararg3:  
%	vararg4:  
%	vararg5:  
%	vararg6:  
%

% ----------------------------- Copyright (C) 2016 Bob Spunt -----------------------------
%	Created:  2016-05-16
%	Email:     spunt@caltech.edu
% ________________________________________________________________________________________
def = { ...
'rmvoxoutlier',         0               , ...
'rmsuboutlier',         0               , ...
'gridflag',             0               , ...
'savefig_combined',     0               , ...
'savefig_separate',     0               , ...
'grid_size',            []              , ...
'row_labels',           ''              , ...
'col_labels',           ''              , ...
'groupspace',           .5              , ...
'tag',                  ''              , ...
'xlab',                 ''              , ...
'roinames',             ''              , ...
'conpat',               'con*nii*'              , ...
'ylab',                 'Percent Signal Change'             , ...
'figpos',               [.05 .05 .14 .85]               , ...
'axgap_h',              .08             , ...
'axgap_v',              .12            , ...
'margin_lr',            [.10 .05]               , ...
'margin_lowerupper',    [.05 .05]               , ...
'basefontsize',         10              , ...
'fontname',             'Helvetica Neue'  ...
	};
vals = setargs(def, varargin);
if nargin < 4, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end

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
if all([isempty(row_labels) nrow>1]), row_labels = strcat({'Group '}, cellstr(num2str((1:nrow)'))); end
if all([isempty(col_labels) ncol>1]), col_labels = strcat({'Level '}, cellstr(num2str((1:ncol)'))); end
if isempty(roinames), [~,roinames] = cellfileparts(roifiles); end
if isempty(grid_size), grid_size = [nroi 1]; end
if all([length(row_labels)~=nrow nrow>1]), printmsg('# of row_labels must match # of rows in conidx'); return; end
if all([length(col_labels)~=ncol ncol>1]), printmsg('# of col_labels must match # of rows in conidx'); return; end
if length(roinames)~=nroi, printmsg('# of roinames must match # of rois'); return; end
if prod(grid_size)~=nroi, printmsg('Grid size must match # of rois'); return; end

% | GET CONTRAST FILES
% | ======================================================================================
conpat = fullfile(level1dirs, conpat);
allcon = cellfun(@files, conpat, 'Unif', false); 
allcon = cellfun(@(x) x(conidx(:)), allcon, 'Unif', false);
allcon = vertcat(allcon{:});
if ncond < max(conidx(:))
    tmpidx = zeros(ncond, 2);
    tmpidx(:,1) = 1:ncond; 
    tmpidx(:,2) = conidx(:); 
    tmpidx = sortrows(tmpidx, 2);
    conidx(tmpidx(:,1)) = 1:ncond; 
end
conidx(:) = 1:ncond; 

% | READ IMAGE VOLUME FILES
% | ======================================================================================
[CON, CONHDR] = bspm_read_vol(allcon, 'reshape');
[ROI, ROIHDR] = bspm_read_vol(roifiles, 'reshape');
for i = 1:nroi
    tmpdata      = CON(ROI(:,i)>0,:);
    if rmvoxoutlier, tmpdata  = bob_outlier2nan(tmpdata, rmvoxoutlier); end
    roidata{i}   = reshape(nanmean(tmpdata), ncond, nlevel1)';
    if rmsuboutlier, roidata{i} = bob_outlier2nan(roidata{i}, rmsuboutlier); end
end

% | UMMMMMMM.... PLOT!
% | ======================================================================================
if gridflag
    hfig         = figure('Color', 'white', 'Units', 'normal', 'Position', figpos);
    [hpan, pidx] = pgrid(grid_size(1), grid_size(2), 'parent', hfig, 'backg', 'white'); 
end
for i = 1:nroi
 
    cg      = row_labels; 
    legname = col_labels;
    cyl     = ylab;
    cxl     = xlab; 
    if gridflag
        cyl = '';
        if i < nroi, cxl = ''; end
        if i==1, legname = ''; end
    end

    % | Plot
    if gridflag, axes('parent', hpan(i)); end

    if ~gridflag, hfig(i) = figure('Color', 'white', 'Units', 'normal', 'Position', figpos); end
    
    hall{i} = barpatch(roidata{i}, ...
            'newfig',   0, ...
            'xl',       cxl, ...
            'yl',       cyl, ...
            'groupidx', conidx, ...
            'fontname', fontname, ...
            'barname', legname, ...
            'fontsize', basefontsize, ...
            'groupname', cg, ...
            'groupspace', groupspace, ...
            't', roinames{i}, ...
            'yticklength', 5);

end

% | CLEANUP
% | ======================================================================================
try
    ht = findall(hfig, 'tag', 'title');
    set(ht, 'fontsize', get(ht(1), 'fontsize')*1.25);
catch
end
try
    hl = findall(hfig, 'tag', 'legend');
    set(hl, 'fontsize', get(hl(1), 'fontsize')*1.15);
catch
end

% | SAVE
% | ======================================================================================
if savefig_combined
    outname = sprintf('ROIPLOT_%s_%dSDs.png', tag, rmsuboutlier);
    fprintf('\n\nexport_fig(''%s'')\n\n', outname); 
end

if savefig_separate
    for i = 1:nroi
        dvlab = sprintf('%s-%s', tag, regexprep(roinames{i}, '\s', '_'));
        outstr = {sprintf('ROIPLOT_%s_', dvlab) sprintf('SubOut-%dSD_', rmsuboutlier) sprintf('VoxOut-%dSD', rmvoxoutlier)};
        outname = strcat(outstr{find([1 rmsuboutlier rmvoxoutlier])});
        outname = regexprep(outname, '_$', '');
        export_fig(strcat(outname, '.jpg'), '-q100', hfig(i));
        cmd{i} = sprintf('export_fig(strcat(%s, ''.jpg''), ''-q100'', %s', outname, handle2str(hfig(i))); 
    end
    disp(cmd)
end

% | OUTPUT
% | ======================================================================================
if nargout==1
    if gridflag
        varargout{1} = hpan; 
    else
        varargout{1} = hfig;
    end
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
