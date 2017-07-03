function matlabbatch = bspm_level2_tstt_batch(condirpat1, condirpat2, varargin)
% BSPM_LEVEL2_TSTT_BATCH
%
%   USAGE: bspm_level2_tstt_batch(analpat1, analpat2, labels, outaffix, mask, conidx)
%
%   ARGUMENTS:
%

% ----------------------------------- Copyright (C) 2014 -----------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
def =   { ...
        'conpat',       'con*.nii*',                ...
        'conidx',       [],                     ...
        'grouplabels',  {'Group1' 'Group2'},    ...
        'tag',          [],                     ...
        'implicitmask',  0,                     ...
        'mask',         bspm_greymask,          ...
        'omitpat',      {''},   ...
        'pctgroup',     90,                     ...
        'runit',        1,                      ...
        'negativecon',   0,                     ...
        'nan2zero',      1,                     ...
        'covariates',   []                      ...
        };
vals = setargs(def, varargin);
if nargin<2, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals);
[cpath, conlist] = files(fullfile(condirpat1, conpat));
if isempty(cpath), error('No contrasts found!'); end
[ucon, idx2ex, idx2sub] = unique(conlist);
if ~isempty(conidx), idx2ex = idx2ex(conidx); end
conname = bspm_con2name(cpath(idx2ex));
conpats = conlist(idx2ex);
if ~runit, allinput = cell(size(conname)); end
for c = 1:length(conname)
    con1 = files(fullfile(condirpat1, conpats{c}));
    con2 = files(fullfile(condirpat2, conpats{c}));
    fprintf('\n| Working on %s, %d of %d', conname{c}, c, length(conname));
    if runit
        bspm_level2_tstt(con1, con2, 'omitpat', omitpat, 'grouplabels', grouplabels, 'tag', tag, 'implicitmask', implicitmask, 'mask', mask, 'pctgroup', pctgroup, 'negativecon', negativecon, 'nan2zero', nan2zero, 'covariates', covariates);
    else
        allinput{c} = bspm_level2_tstt(con1, con2, 'omitpat', omitpat, 'grouplabels', grouplabels, 'tag', tag, 'implicitmask', implicitmask, 'mask', mask, 'pctgroup', pctgroup, 'negativecon', negativecon, 'nan2zero', nan2zero, 'covariates', covariates);
    end
end
fprintf('\n - DONE - \n\n');
end








