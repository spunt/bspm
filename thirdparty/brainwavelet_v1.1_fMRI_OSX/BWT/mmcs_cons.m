function Nbmat = mmcs_cons(M,padM,wlh,whh,bin)
%
% FUNCTION:     mmcs_cons -- Finds chain coefficients (conservative
%                            method).
%                            
% USAGE:        Nbmat = mmcs_mod(M,padM,wlh,whh,bin)
%
% Inputs:       M        -- Matrix to be searched for chains.
%                           This should be a binary matrix indicating
%                           maxima and minima.
%               padM     -- Padded matrix to accommodate neighbourhood
%                           search. E.g. if matrix dimensions are
%                           3x3x3, and neighbourhood for searching is
%                           2x1 either side of coefficients, padM
%                           should have dimensions 7x5x3.
%               wlh      -- Width of neighbourhood to one side of
%                           the coefficient, in y dimension.
%               whh      -- Height of neighbourhood to one side of the
%                           coefficient, in x dimension.
%               bin      -- Binary flag [Default=1]. Tells function
%                           whether to output the number of neighbours
%                           for a given coefficient [0], or a binary
%                           digit indicating the presence/absence of
%                           neighbours [1].
%
% Output:       Nbmat    -- Output neighbourhood matrix.
%
% EXAMPLE:      NBmat = mmcs_mod(3Dmat,3Dpadmat,2,1,0)
%
% AUTHOR:       Ameera X Patel
% CREATED:      02-01-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     7
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: mmcs_mod.m 7 30-01-2014 BWTv1.1 axpatel


%% check inputs (matrix dimensions)

fname=mfilename;
if nargin<1
    help(fname); return
end

[NX,Nscales,Nts]=size(M);
[NXx,Nscalesx,Ntsx]=size(padM);

if ( NXx~=(NX+(2*wlh)) || Nscalesx~=(Nscales+(2*whh)) || Ntsx~=Nts )
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Invalid window assignment or matrix dimensions\n')
    Nbmat=[];
    return
end

%% check inputs (number,types,values)

ok_whh=1:(Nscalesx-Nscales)/2;
ok_wlh=1:(NXx-NX)/2;

err=struct();
err(1).inp=chkninput(nargin,[4,5],nargout,[0,1],fname);

err(2).inp=chkintype(M,'double',fname);
err(3).inp=chkintype(padM,'double',fname);
err(4).inp=chkintype(wlh,'numeric',fname,ok_wlh);
err(5).inp=chkintype(whh,'numeric',fname,ok_whh);

if exist('bin','var')
    err(6).inp=chkintype(bin,'numeric',fname,{'0','1'});
else
    err(6).inp=0;
    bin=1;
end

if sum(cat(1,err.inp))>=1
    Nbmat=[]; return
end

%% binarize input matrix M

M(M~=0)=1;

%% neighbourhood search scales 1,2

Nbmat=zeros(NX,Nscales,Nts);
Nbmat(:,1,:)=nghbrcount(padM(:,(1:1+whh+whh),:),wlh,whh,NX,1,Nts,bin);

binNbmat=Nbmat(:,1,:);
binNbmat(binNbmat>0)=1;

mmcs_coef=zeros(NX,1,Nts);
mmcs_coef((binNbmat+M(:,1,:))==2)=1;

mmcmat=padM;
mmcmat(:,1+whh,:)=circularize(mmcs_coef,'y',wlh);

mmcmat(:,(3+whh:Nscalesx),:)=0;

for i=2:Nscales
    mmcmat(:,i+whh,:)=padM(:,i+whh,:);
    Nbmat(:,i,:)=nghbrcount(mmcmat(:,(i:i+whh+whh),:),wlh,whh,NX,1,Nts,bin);
    binNbmat=Nbmat(:,i,:);
    binNbmat(binNbmat>0)=1;
    mmcs_coef=zeros(NX,1,Nts);
    mmcs_coef((binNbmat+M(:,i,:))==2)=1;
    mmcmat(:,i+whh,:)=circularize(mmcs_coef,'y',wlh);
    
end

end