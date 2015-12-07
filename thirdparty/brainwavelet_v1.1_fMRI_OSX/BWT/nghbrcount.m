function Nbmat = nghbrcount(padM,wlh,whh,NX,Nscales,Nts,bin)
%
% FUNCTION:     nghbrcount -- Identifies positive neighbours within a
%                             given neighbourhood.
%                            
% USAGE:        Nbmat = nghbrcount(padM, wlh, whh, NX, Nscales, Nts, bin)
%
% Inputs:       padM       -- Padded matrix to accommodate neighbour
%                             search. E.g. if matrix dimensions are
%                             3x3x3, and neighbourhood for searching is
%                             2x1 either side of coefficients, padM
%                             should have dimensions 7x5x3.
%               wlh        -- Width of neighbourhood to one side of
%                             the coefficient, in y dimension.
%               whh        -- Height of neighbourhood to one side of
%                             the coefficient, in x dimension.
%               NX         -- Size of original matrix in y direction.
%               Nscales    -- Size of original matrix in x direction.
%               Nts        -- Size of original matrix in z direction.
%               bin        -- Binary flag [Default=1]. Tells function
%                             whether to output the number of neighbours
%                             for a given coefficient [0], or a binary
%                             digit indicating the presence/absence of
%                             neighbours [1].
%
% Output:       Nbmat      -- Output neighbourhood matrix
%
% EXAMPLE:      Nbmat = nghbrcount(3Dpadmat,2,1,100,6,1000,1)
%
% AUTHOR:       Ameera X Patel
% CREATED:      31-12-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     6
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: nghbrcount.m 6 30-01-2014 BWTv1.1 axpatel


%% check inputs (matrix dimensions)

fname=mfilename;
if nargin<1
    help(fname); return
end

[NXx,Nscalesx,Ntsx]=size(padM);

if ( NXx~=(NX+(2*wlh)) || Nscalesx~=(Nscales+(2*whh)) || Ntsx~=Nts )
    cprintf('err','\nInvalid window assignment or matrix dimensions.\n')
    Nbmat=[];
    return
end

%% check inputs (number,types,values)

ok_whh=1:(Nscalesx-Nscales)/2;
ok_wlh=1:(NXx-NX)/2;

err=struct();
err(1).inp=chkninput(nargin,[6,7],nargout,[0,1],fname);

err(2).inp=chkintype(padM,'double',fname);
err(3).inp=chkintype(wlh,'numeric',fname,ok_wlh);
err(4).inp=chkintype(whh,'numeric',fname,ok_whh);
err(5).inp=chkintype(NX,'numeric',fname);
err(6).inp=chkintype(Nscales,'numeric',fname);
err(7).inp=chkintype(Nts,'numeric',fname);

if exist('bin','var')
    err(8).inp=chkintype(bin,'numeric',fname,{'0','1'});
else
    err(8).inp=0;
    bin=1;
end

if sum(cat(1,err.inp))>=1
    Nbmat=[]; return
end

%% binarize input matrix

padM(padM~=0)=1;

%% generate window-shifted matrices and sum

% whh=(height-1)/2;
% wlh=(width-1)/2;

ctx=1;
outmat=zeros(NX,Nscales,Nts);

for i=-wlh:wlh
    for ii=-whh:whh
        eval(['Y' num2str(ctx) ...
            '=padM(1+wlh+i:NX+wlh+i,1+whh+ii:Nscales+whh+ii,:);']);
        outmat=outmat+eval(['Y' num2str(ctx)]);
        ctx=ctx+1;
        clear Y*;
        
    end
end

outmat=outmat-padM(1+wlh:NX+wlh,1+whh:Nscales+whh,:);

% list=who('Y*');
% list(strcmp(list,strcat('Y',num2str(ceil(2*whh*wlh+wlh+whh+0.5)))))=[];

%% output neighbour matrix

if bin==0;
    Nbmat=outmat;
elseif bin==1
    Nbmat=(outmat>0)+0;
end

end