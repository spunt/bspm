function [] = WriteOutNii(ts,outfname,Info)
%
% FUNCTION:     WriteOutNii -- writes NIfTI file of input compressed or
%                              uncompressed 2D matrix given dimensions and
%                              compression data.
%                            
% USAGE:        WriteOutNii(ts,outfname,Info)
%
% Inputs:       ts          -- 2D matrix of time series to write out as 
%                              NIfTI.
%               outfname    -- Output file prefix.
%               Info        -- Structure containing the following fields:
%                               hdr   - Containing dataset header.
%                               tsInd - Indices of time series in original
%                                       uncompressed 2D matrix].
%                               dim   - Original dataset dimensions.
%
%                               If dataset has been compressed, tsInd and
%                               dim are required fields.
%
% Output:       This function will output matrix ts as a NIfTI file, with 
%               name outfname.nii.gz, into the current directory.
%
% EXAMPLE:      WriteOutNii(2Dmat,'MyNIfTI.nii.gz',Info)
%
% AUTHOR:       Ameera X Patel
% CREATED:      14-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     2
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: WriteOutNii.m 2 30-01-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return;
end

if chkninput(nargin,[3,3],nargout,[0,0],fname)>=1
    return;
end

err=struct();
err(1).inp=chkintype(ts,'double',fname);
err(2).inp=chkintype(outfname,'char',fname);
err(3).inp=chkintype(Info,'struct',fname);

if sum(cat(1,err.inp))>=1;
    return
end

%% check hdr field exists

if (isfield(Info,'hdr')&&isstruct(Info.hdr))==0
    cprintf('_err','*** BrainWavelet Error ***')
    cprintf('err','cannot locate dataset hdr field');
end

%% Uncompress matrix / orient ts

tsd=size(ts);
hdrDim=Info.hdr.dime.dim;

if isfield(Info,'tsInd')
    CompL=length(Info.tsInd);
    if not(isfield(Info,'dim'))
       cprintf('_err','*** BrainWavelet Error ***\n')
       cprintf('err','Cannot ascertain dimensions of original dataset.\n');
       cprintf('err','Info field indicates the dataset was compressed\n');
       cprintf('err','for speed. Cannot uncompress without dimensions.\n');
       return;
    end
    if CompL==tsd(2);
       ts=ts';
       tsd=size(ts);
    elseif CompL~=tsd(1) || length(tsd)>2
        cprintf('_err','*** BrainWavelet Error ***\n')
        cprintf('err','hdr/dataset dimension mismatch.\n');
        return;
    end
    UnCts=zeros(Info.dim(1)*Info.dim(2)*Info.dim(3),tsd(2));
    UnCts(Info.tsInd,:)=ts;
    tsd=size(UnCts);
end

%% check ts-hdr nvox dimension match in uncompressed matrix, and orient ts

nvox=hdrDim(2)*hdrDim(3)*hdrDim(4);
if nvox==tsd(2);
	ts=ts';
	tsd=size(ts);
	UnCts=ts;
elseif nvox~=tsd(1) || length(tsd)>2
	cprintf('_err','*** BrainWavelet Error ***\n')
	cprintf('err','hdr/dataset dimension mismatch.\n');
	return;
end

%% reshape to 4D

tsRs=reshape(UnCts,hdrDim(2),hdrDim(3),hdrDim(4),tsd(2));

%% edit hdr for NBriks

if tsd(1)==1;
    Info.hdr.dime.dim(1)=3;
    Info.hdr.dime.dim(5)=1;
elseif tsd(2)~=Info.hdr.dime.dim(5);
    Info.hdr.dime.dim(5)=tsd(2);
end

%% Write NIfTI

OutNii=struct();
OutNii.hdr=Info.hdr;
OutNii.img=tsRs;

save_nii(OutNii,outfname);

end