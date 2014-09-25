function bspm_omit_vols(folderpat, omitpat)
% BSPM_OMIT_VOLS
%
% USAGE: bspm_omit_vols(folderpat, omitpat)  
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
%       >> omitpat = {'fad*1-000001*nii' 'fad*2-000002*nii'};
%       >> bspm_omit_vols(folderpat, omitpat)
%

% ---------------------- Copyright (C) 2014 ----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, error('USAGE: bspm_omit_vols(folderpat, omitpat)'); end

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
 
 
 
 
