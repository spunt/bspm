function [ts,Info,error] = ParseInNii(InNii,varargin)
%
% FUNCTION:     ParseInNii -- Reads NIfTI file and converts it to a 
%                             compressed or uncompressed 2D matrix.
%                            
% USAGE:        ParseInNii(InNii,varargin)
%
% Inputs:       InNii      -- Input file name. Input must be a string
%                             containing path to file if required.
%
%               Additional Input Options:
%               (These must be specified as MATLAB string-value pairs).
%
%               compress   -- Binary flag indicating whether to compress 
%                             the dataset for speed [1] or not [0]. 
%                             [Default=1].
%
% Outputs:      ts         -- 2D matrix representing time x voxels.
%               Info       -- Structure containing the following fields:
%                               hdr   - Containing dataset header.
%                               tsInd - Indices of time series in original
%                                       uncompressed 2D matrix].
%                               dim   - Original dataset dimensions.
%
%                               If dataset has been compressed, tsInd and
%                               dim will be automatically output.
%
% EXAMPLE:      [ts,Info] = ParseInNii('MyNIfTI.nii.gz','compress',1)
%
% AUTHOR:       Ameera X Patel
% CREATED:      14-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     8
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: ParseInNii.m 8 10-02-2014 BWTv1.1 axpatel

%% parse options and check files/inputs

error=0;
fname=mfilename;
if nargin<1
    help(fname); return;
end

if chkninput(nargin,[1,3],nargout,[0,3],fname)>=1
    ts=[]; if nargout>=2; Info=[]; end
    return;
end

DefaultOpts=struct('compress',1);
Opts=parseInOpts(DefaultOpts,varargin);

if exist(InNii,'file')~=2;
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Cannot find input NIfTI file.\n\n')
    error=error+1; ts=[]; if nargout>=2; Info=[]; end
    return;
end

err=struct();
err(1).inp=chkintype(InNii,'char',fname);
err(2).inp=chkintype(Opts.compress,'numeric',fname,{'0','1'});

if sum(cat(1,err.inp))>=1;
    error=error+1; ts=[]; if nargout>=2; Info=[]; end
    return;
end

%% load matrix into memory

[nii,obQ]=load_nii(InNii);
if obQ==1;
    nii=load_untouch_nii(InNii);
end
M=nii.img;

%% reshape (and compress) matrix

Info=struct();
Info.hdr=nii.hdr;

[x,y,z,t]=size(M);
ts=reshape(M,x*y*z,t);
Info.dim=[x y z t];

if Opts.compress==1
    zInd=all(ts==0,2);
    Info.tsInd=find(zInd==0);
    ts(zInd,:)=[];
end

%% flip matrix dimension for wavelets

ts=double(ts)';

end