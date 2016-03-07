function [D altOut] = extractTS(fn, seeds, opt)
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
%%%         {'PCC'  [  0 -53  26] [8]} ...
%%%         {'mPFC' [  0  52 -06] [8]} ...
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
% keyboard;
if nargin == 2
    opt = 1;
end

if isstruct(fn)
    V1 = fn;
else
    V1 = spm_vol(fn);
end

[tmp XYZmm] = spm_read_vols(V1(1));

if nargin>2  
    load([fileparts(fn) '/' opt]);
    [row, col] = find(R==1);
    ind = setdiff(1:numel(V1),row);
   
    M2 = FastRead(V1(ind));
else
     M2 = FastRead(V1);
end

for zz = 1:length(seeds)
    if isnumeric(seeds{zz}{2});
        rad = seeds{zz}{3}/2;
        dist = sqrt(sum((XYZmm' - repmat(seeds{zz}{2},size(XYZmm,2),1)).^2,2));
        ind = find(dist <= rad);
        if numel(ind)==0;
            ind = find(dist==min(dist));
        end
        D.(seeds{zz}{1}).fn = ['Mask_' num2str(seeds{zz}{3}) 'mm_' num2str(seeds{zz}{2}(1)) '_' num2str(seeds{zz}{2}(2)) '_' num2str(seeds{zz}{2}(3)) '.nii'];
    elseif ischar(seeds{zz}{2});
        h = spm_vol(seeds{zz}{2});
        mi = resizeVol(h,V1(1));
        ind = find(mi==1);
        [aa1 aa2] = fileparts(seeds{zz}{2});
        D.(seeds{zz}{1}).fn = ['Mask_' aa2];
    else
        warning(['Seed #' num2str(ii) ' was not specified correctly.']);
        continue
    end
    
    dd = XYZmm(:,ind)';
    clear dd2;
    [dd2(:,1) dd2(:,2) dd2(:,3)] = ind2sub(V1(1).dim,ind);
    
    D.(seeds{zz}{1}).mni_loc = dd;
    D.(seeds{zz}{1}).vox_loc = dd2;
    D.(seeds{zz}{1}).vec_loc = ind;
    D.(seeds{zz}{1}).ts = [];
    
    %if opt == 1
        y = M2(:,ind);
        D.(seeds{zz}{1}).tsm = nanmean(y,2);
        D.(seeds{zz}{1}).allVox = y;
        D.(seeds{zz}{1}).tsv = QuickSVD(y);
        
        altOut{1}{zz} = seeds{zz}{1};
        altOut{2}(:,zz) = nanmean(y,2);
    %end
    
end
