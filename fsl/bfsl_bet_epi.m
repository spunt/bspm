function bfsl_bet_epi(epi, f, outprefix)
% BFSL_BET_EPI
%
% USAGE: bfsl_bet_epi(epi, f, outprefix)
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
if nargin<3, outprefix = 'b'; end
if nargin<2, f = 0.3; end
if nargin<1, mfile_showhelp; return; end
if ischar(epi), epi = cellstr(epi); end

% | skull strip first image
fprintf('\n | - Running BET on First Image'); 
input   = epi{1};
[p n e] = fileparts(input);
output  = [p filesep outprefix n '.nii.gz'];
command = sprintf('bet %s %s -f %2.1f',input, output, f);
system(command);
gunzip(output);
delete(output);

% | get indices of bad voxels
output  = [p filesep outprefix n '.nii'];
mout    = bspm_read_vol(output, 'reshape'); 
mout    = mout < (nanmean(mout)/8);

% | get all vols
fprintf('\n | - Reading All Images'); 
[d,h]       = bspm_read_vol(epi, 'reshape');
hout        = h;
d(mout, :)  = NaN;
[p,n,e]     = cellfun(@fileparts, epi, 'Unif', false);
outname     = strcat(p, filesep, outprefix, n, e); 

% | loop over volumes
fprintf('\n | - Writing Images: '); 
for i = 1:length(epi)
    printcount(i, length(epi)); 
    hout(i).fname = outname{i}; 
    spm_write_vol(hout(i), reshape(d(:,i), hout(i).dim)); 
end
fprintf('\n | - Done!\n');  
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
