function bspm_bet_epi_batch(epi,f,outprefix,waitbar)
% BSPM_BET_EPI_BATCH
%
% USAGE: bspm_bet_epi_batch(epi,f,outprefix)
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

if nargin<4, waitbar = 0; end
if nargin<3, outprefix = 'b'; end
if nargin<2, f = 0.3; end
if ischar(epi), epi = cellstr(epi); end

% skull strip first image
input = epi{1};
[p n e] = fileparts(input);
output = [p filesep outprefix n '.nii.gz'];
command = sprintf('bet %s %s -f %2.1f',input, output, f);
system(command);
gunzip(output);
delete(output);

% get indices of bad voxels
output = [p filesep outprefix n '.nii'];
v = spm_vol(output);
d = spm_read_vols(v);
badidx = d < (mean2(d)/8);

% loop over remaining volumes
nim = length(epi);
if waitbar
h = waitbar(0,sprintf('Completed 1 of %d Volumes',nim));
end
for i = 1:nim
    input = epi{i};
    [p n e] = fileparts(input);
    output = [p filesep outprefix n '.nii'];
    v = spm_vol(input); 
    d = spm_read_vols(v);
    d(badidx) = NaN;
    v.fname = output; 
    spm_write_vol(v,d);
    if waitbar
    waitbar(i/nim,h,sprintf('Completed %d of %d Volumes',i,nim));
    end
end
if waitbar
close(h);
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
