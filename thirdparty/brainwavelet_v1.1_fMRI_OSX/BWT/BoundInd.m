function Ind = BoundInd(Lj,NX,NJ)
%
% FUNCTION:     BoundInd -- Finds boundary coefficient start/end indices.
%                            
% USAGE:        Ind = BoundInd(Lj,NX,NJ)
%
% Inputs:       Lj       -- Vector of filter widths at various scales.
%               NX       -- Number of time points in time series.
%               NJ       -- Number of scales.
%
% Output:       Ind      -- Vector (dimensions 2 x NJ) that contains start
%                           and end indices of boundary coefficients.
%
% EXAMPLE:      indices = BoundInd(Lj,100,7)
%
% AUTHOR:       Ameera X Patel
% CREATED:      11-01-2014
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     5
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: BoundInd.m 5 30-01-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); 
    Ind=[]; return;
end

if chkninput(nargin,[3,3],nargout,[0,1],fname) >=1
    Ind=[]; return;
end

err=struct();
err(1).inp=chkintype(Lj,'double',fname);
err(2).inp=chkintype(NX,'numeric',fname);
err(3).inp=chkintype(NJ,'numeric',fname);

if sum(cat(1,err.inp))>=1;
    Ind=[]; return;
end

clear fname err

Ljd=size(Lj);

if sum(Ljd>=NJ)==0
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Not enough filter info for number of scales\n');
    cprintf('err','Please ensure that Vj is at least of length %s\n',...
        num2str(NJ));
    Ind=[]; return;
end

if isvector(Lj)==0
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Please ensure that Lj is a vector\n');
    Ind=[]; return;
end

%% find boundary indices

Ind=zeros(2,NJ);
Ind(2,:)=Lj(1:NJ)-2;
Ind=Ind+1;
Ind(Ind>NX)=NX;

end