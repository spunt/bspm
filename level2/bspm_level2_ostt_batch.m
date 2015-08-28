function allinput = bspm_level2_ostt_batch(condirpat, varargin)
% BSPM_LEVEL2_OSTT_BATCH
%
%   USAGE: allinput = bspm_level2_ostt_batch(condirpat, varargin)
%
%   ARGUMENTS:
%       condirpat: pattern for finding level 1 dirs containing con images
%

% ------------------------------------- Copyright (C) 2014 -------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
def = { 'outdir',       [],     ...
        'conidx',       [],     ...     
        'tag',          [],     ...
        'implicit',     0,      ...
        'mask',         '',     ...
        'pctgroup',     90,     ...
        'negativecon',  0,      ...
        'conpat',       'con*nii', ...
        'runit',        1,      ...
        'nan2zero',     1       ...
        };
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
fprintf('\n\t= CURRENT SETTINGS =\n'); disp(vals); 
if ~regexp(condirpat, pwd), condirpat = fullfile(pwd, condirpat); end
[cpath, conlist] = files(fullfile(condirpat, conpat)); 
if isempty(cpath), error('No contrasts found!'); end
[ucon, idx2ex, idx2sub] = unique(conlist);
if ~isempty(conidx), idx2ex = idx2ex(conidx); end
conname = bspm_con2name(cpath(idx2ex)); 
for c = 1:length(conname)
    cons = cpath(idx2sub==idx2ex(c)); 
    fprintf('\n| Working on %s, %d of %d', conname{c}, c, length(conname));
    if runit
        bspm_level2_ostt(cons, 'outdir', outdir, 'tag', tag, 'implicit', implicit, 'mask', mask, 'pctgroup', pctgroup, 'negativecon', negativecon, 'nan2zero', nan2zero); 
    else
        allinput{c} = bspm_level2_ostt(cons, 'outdir', outdir, 'tag', tag, 'implicit', implicit, 'mask', mask, 'pctgroup', pctgroup, 'negativecon', negativecon, 'nan2zero', nan2zero); 
    end
end
fprintf('\n - DONE - \n\n'); 


 
 
 
 
