function [] = bspm_level2_ostt_batch(condirpat, minN, mask, rmoutliers, conidx, omit, conpat)
% BSPM_LEVEL2_OSTT_BATCH
%
%   USAGE: bspm_level2_ostt_batch(condirpat, minN, mask, rmoutliers, conidx, omit, conpat)
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

if nargin<7, conpat='con*img'; end
if nargin<6, omit= []; end
if nargin<5, conidx= []; end
if nargin<4, rmoutliers = 0; end
if nargin<3, mask = []; end
if nargin<2, minN=[]; end
if nargin<1, disp('USAGE: bspm_level2_ostt_batch(condirpat, minN, mask, rmoutliers, conidx, omit, conpat)'); return; end
condirpat = [pwd filesep condirpat];
condirs = files(condirpat);
if isempty(condirs), bob_display_message('No directories found! Check working dir and contrast dir pattern'); return; end
[cpath, conlist] = files([condirs{1} filesep conpat]);
if conidx, conlist = conlist(conidx); end
for c = 1:length(conlist)
    cons = files([condirpat filesep conlist{c}]);
    bob_display_message(sprintf('Working on Contrast %d of %d: %s', c, length(conlist), conlist{c}));
    if ~isempty(omit)
        bspm_level2_ostt(cons, minN, mask, rmoutliers, omit);
    else 
        bspm_level2_ostt(cons, minN, mask, rmoutliers);
    end
        
end
    


 
 
 
 
