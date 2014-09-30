function [] = bspm_level2_ostt_batch_spm(condirpat, name, mask, implicitTAG, conidx, conpat)
% BSPM_LEVEL2_OSTT_BATCH
%
%   USAGE: bspm_level2_ostt_batch(condirpat, name, mask, implicitTAG, conidx, conpat)
%
%   ARGUMENTS:
%       condirpat: pattern for finding level 1 dirs containing con images
%       name: affix for directory in _groupstats_/<analysis_name> 
%       (optional) mask: mask file to use (default = bspm_greymask)
%       (optional) implicitTAG: 0 = no implicit masking; 1 = yes (default)
%       (optional) conidx: indices for contrasts to include
%       (optional) conpat: indices for pattern for finding contrasts
%

% ------------------------------------- Copyright (C) 2014 -------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<6, conpat='con*img'; end
if nargin<5, conidx=[]; end
if nargin<4, implicitTAG=1; end
if nargin<3, mask = bspm_greymask; end
if nargin<2, disp('USAGE: bspm_level2_ostt_batch(analysispat, name, mask, implicitTAG, conidx, conpat)'); return; end
condirs = files(condirpat);
if isempty(condirs), bob_display_message('No directories found! Check working dir and contrast dir pattern'); return; end
[cpath,conlist] = files([condirs{1} filesep conpat]);
if conidx, conlist = conlist(conidx); end
for c = 1:length(conlist)
    cons = files([condirpat filesep conlist{c}]);
    bob_display_message(sprintf('Working on Contrast %d of %d: %s', c, length(conlist), conlist{c}));
    bspm_level2_ostt_spm(cons, name, mask, implicitTAG);
end
    


 
 
 
 
