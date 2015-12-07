function inds = ParseInds(NTS,NBatch)
%
% FUNCTION:     ParseInds -- Computes indices of voxels to process in
%                            batch mode.
%                            
% USAGE:        ParseInds(NVox,NBatch)
%
% Inputs:       NTS       -- Number of time series.
%               NBatch    -- Number of batches in which to split the time
%                            series.
%
% Output:       inds      -- Start and end indices of batches. Output 
%                            matrix has dimensions 2 x NBatch.
%
% EXAMPLE:      inds = ParseInds(32155,3)
%
% AUTHOR:       Ameera X Patel
% CREATED:      22-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     4
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: ParseInds.m 4 30-01-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return;
end

err=struct();
err(1).inp=chkninput(nargin,[2,2],nargout,[0,1],fname);
err(2).inp=chkintype(NTS,'numeric',fname);
err(3).inp=chkintype(NBatch,'numeric',fname);

if sum(cat(1,err.inp))>=1
    inds=[]; return;
end

%% find indices

step=round(NTS/NBatch);

inds=zeros(2,NBatch);
inds(2,:)=(1:NBatch).*step;
inds(2,NBatch)=NTS;
inds(1,:)=inds(2,:)-step+1;

end