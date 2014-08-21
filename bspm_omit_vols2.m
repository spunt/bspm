function bspm_omit_vols2(input)
% BSPM_OMIT_VOLS
%
% USAGE: bspm_omit_vols(folderpat, omitpat, fidx)  
%
%   INPUTS:
%       folderpat:  for defining the group of folders
%       omitpat:    for defining files within each folder
%
%   EXAMPLE:
%       This omits the volumes matching pattern '*nii' from
%       all folders matching the pattern '/data/my_study/sub*/epi'
%
%       >> folderpat = '/data/my_study/sub*/epi';
%       >> omitpat = '*nii';
%       >> fidx = 1:4;
%       >> bspm_omit_vols(folderpat, omitpat)
%

% ---------------------- Copyright (C) 2014 ----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014



folderpat = input.folderpat;
omitpat = input.omitpat;

% get folders
foldernames = files(folderpat);
nfold = length(foldernames);

% move files
for i = 1:nfold
    outdir = [foldernames{i} filesep '_omit_'];
    if ~isdir(outdir), mkdir(outdir); end
    for f = 1:length(omitpat)
        p = files([foldernames{i} filesep omitpat{f}]);
        if ~isempty(p)
            movefile(char(p),outdir);
        end
    end
end
 
 
 
 
