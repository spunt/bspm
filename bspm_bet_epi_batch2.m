function bspm_bet_epi_batch2(input)
% BSPM_BET_EPI_BATCH2
%
% USAGE: bspm_bet_epi_batch2(input.epipat)
%

% ------------ Copyright (C) 2014 ------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

outprefix = 'b';
f = 0.3;
epipat = input.epipat;
if iscell(epipat), epipat = char(epipat); end
epi = files(epipat);

% skull strip first image
input = epi{1};
[p n e] = fileparts(input);
output = [p filesep outprefix n '.nii.gz'];
command = sprintf('system(''bet %s %s -f %2.2f'')',input, output, f);
eval(command);
gunzip(output);
delete(output);

% get indices of bad voxels
output = [p filesep outprefix n '.nii'];
v = spm_vol(output);
d = spm_read_vols(v);
badidx = d < (mean2(d)/8);

% loop over remaining volumes
input = [];
nim = length(epi);
for i = 1:nim
    input = epi{i};
    [p n e] = fileparts(input);
    output = [p filesep outprefix n '.nii'];
    v = spm_vol(input); 
    d = spm_read_vols(v);
    d(badidx) = NaN;
    v.fname = output; 
    spm_write_vol(v,d);
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
