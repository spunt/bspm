function bfsl_bet_epi2(epi, outprefix)
% BFSL_BET_EPI
%
% USAGE: bfsl_bet_epi2(epi, outprefix)
%
% ARGUMENTS
%   
%

% ----------------- Copyright (C) 2014 -----------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<2, outprefix = 'b'; end
if nargin<1, disp('USAGE: bfsl_bet_epi2(epi, outprefix)'); return; end
if iscell(epi), epi = char(epi); end
fprintf('\n | - Running BET'); 
[p, n, e] = fileparts(epi);
output  = [p filesep outprefix n '.nii.gz'];
command = sprintf('bet %s %s -F', epi, output);
system(command);
fprintf('\n | - Done!\n');  
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
