function [] = bspm_level2_tstt_batch(analpat1, analpat2, labels, outaffix, mask, conidx)
% BSPM_LEVEL2_TSTT_BATCH
%
%   USAGE: bspm_level2_tstt_batch(analpat1, analpat2, labels, outaffix, mask, conidx)
%
%   ARGUMENTS:
%       analpat1: pattern for finding analysis directories for group 1
%       analpat2: pattern for finding analysis directories for group 2
%       outaffix: affix for directory in _groupstats_/<analysis_name> 
%       (optional) mask: mask file to use (default = bspm_greymask)
%       (optional) conidx: indices for contrasts to include
%

% ----------------------------------- Copyright (C) 2014 -----------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<6, conidx=[]; end
if nargin<5, mask = bspm_greymask; end
if nargin<4, outaffix = 'TSST'; end
if nargin<3, disp('USAGE: bspm_level2_tstt_batch(analpat1, analpat2, labels, outaffix, mask, conidx)'); return; end
tmpdirs = files(analpat1);
conlist = files([tmpdirs{1} filesep 'con*img'], 'filename');
conlist(conidx) = [];
for c = 1:length(conlist)
    
    cons1 = files([analpat1 filesep conlist{c}]);
    cons2 = files([analpat2 filesep conlist{c}]);
    bspm_display_message(sprintf('Working on Contrast %d of %d: %s', c, length(conlist), conlist{c}));
    bspm_display_message(sprintf('%s: %d   %s:%d', labels{1}, length(cons1), labels{2}, length(cons2)), '-');
    bspm_level2_tstt(cons1, cons2, labels, outaffix, mask);
    
end
    
 
 
 
 
 
 
 
 
