function [D altOut] = extractTS2(fn, seeds,opt)
%%%  Extract a time series from a 4D Nifti file.
%%%
%%% INPUTS:
%%% fn:     Path to 4D nifti file.
%%% seeds:  Cell Array of Inputs { {'Label1' [mni coords] [seed diameter]} ...}
%%% opt:    If opt = 1, the time series will be extracted.
%%%         If opt = 0, only the matrix indcies of the seed will be
%%%         computed.
%%%
%%% OUTPUTS
%%% D.fn =  A new filename based on the location and size of the seed.
%%% D.vox_loc = location of the seed in 3D matrix coordinates.
%%% D.vec_loc = location of the seed in vector coordinates.
%%% D.ts =  The time series as computed by the first eigen vector method.
%%% D.tsm = The time series as computed by a simple mean of the values in
%%% the seed at each time point.
%%%
%%% Example: 
%%%    extractTS(rsscan{1}, { ...
%%%         {'PCC'  [  0 -53  26] [10]} ...
%%%         {'mPFC' [  0  52 -06] [10]} ...
%%%         {'lLPC' [-48 -62  36] [8]} ...
%%%         {'LRPC' [ 46 -62  32] [8]} ...
%%%         });
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Last Updated Dec. 11 2012;
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.
altOut{1} = {};
altOut{2} = [];


if isstruct(fn)
    V1 = fn;
else
    V1 = spm_vol(fn);
end

if nargin>2  
    load([fileparts(fn) '/' opt]);
    [row, col] = find(R==1);
    ind = setdiff(1:numel(V1),row);
   
    V1 = V1(ind);
end

for zz = 1:length(seeds)
    if isnumeric(seeds{zz}{2});
        [ml vi] = getMatCoord(V1(1),seeds{zz}{2},seeds{zz}{3});
    elseif ischar(seeds{zz}{2});
        h = spm_vol(seeds{zz}{2});
        mi = resizeVol(h,V1(1));
        vi = find(mi==1);
        ml = [];
        [ml(:,1), ml(:,2), ml(:,3)] = ind2sub(V1(1).dim,vi);
    else
        warning(['Seed #' num2str(ii) ' was not specified correctly.']);
        continue
    end
    
    mni = [ml ones(size(ml,1),1)]*V1(1).mat';
    
   
    D.(seeds{zz}{1}).mni_loc = mni;
    D.(seeds{zz}{1}).vox_loc = ml;
    D.(seeds{zz}{1}).vec_loc = vi;
    
    X = ExtractNIFTIdata(V1,ml);
    D.(seeds{zz}{1}).tsm = nanmean(X,2);
    D.(seeds{zz}{1}).allVox = X;
    D.(seeds{zz}{1}).tsv = QuickSVD(X(:,~isnan(mean(X))));
    
    altOut{1}{zz} = seeds{zz}{1};
    altOut{2}(:,zz) = nanmean(X,2);
end
