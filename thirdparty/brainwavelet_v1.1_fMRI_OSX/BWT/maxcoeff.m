function localmm = maxcoeff(m,thresh,mcoef)
%
% FUNCTION:     maxcoeff -- Finds local maxima and minima greater than
%                           threshold. Coefficients > mcoef of the local 
%                           maximum/minimum are also included.
%                            
% USAGE:        localmm = maxcoeff(m,thresh,mcoef)
%
% Inputs:       m        -- 3D matrix of wavelet decomposed time series
%                           NX time points x NJ scales x NTS time series.
%               thesh    -- Threshold to cap coefficients.
%               mcoef    -- Maximum/Minimum leniency coefficient.
%                           [Default=0.5].
%
% Output:       localmm  -- Binary matrix of maxima and minima.
%
% EXAMPLE:      localmm = maxcoeff(3Dmat,10,0.5)
%
% AUTHOR:       Ameera X Patel
% CREATED:      26-12-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     4
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: maxcoeff.m 4 30-01-2014 BWTv1.1 axpatel


%% check inputs

fname=mfilename;
if nargin<1
    help(fname); return
end

err=struct();
err(1).inp=chkninput(nargin,[2,3],nargout,[0,1],fname);
err(2).inp=chkintype(m,'double',fname);
err(3).inp=chkintype(thresh,'numeric',fname);

if exist('mcoef','var')
    err(4).inp=chkintype(mcoef,'numeric',fname);
else
    mcoef=0.5;
    err(4).inp=0;
end

if sum(cat(1,err.inp))>=1;
    localmm=[]; return;
end

%% determine dimensions

[NX,Nscales,Nts] = size(m);
localmax  = zeros(NX,Nscales,Nts);
localmin  = zeros(NX,Nscales,Nts);

%% initiate neighbourhood indices

t   = 1:NX;
tp1 = [t(NX) t(1:(NX-1))];
tp2 = [t(NX-1) t(NX) t(1:(NX-2))];
tm1 = [t(2:NX) t(1)];
tm2 = [t(3:NX) t(1) t(2)];

%% find local maxima and minima

for i = 1:Nscales;
    x=m(:,i,:);
    Nxmax=max([x(tm2(:),1,:),x(tm1(:),1,:),x(t(:),1,:),...
        x(tp1(:),1,:),x(tp2(:),1,:)],[],2);
    Nxmax=x./Nxmax;
    localmax(:,i,:) = ((Nxmax>=mcoef)+0);
    localmax(:,i,:) = localmax(:,i,:) .* (abs(x)>thresh);
    
    Nxmin=min([x(tm2(:),1,:),x(tm1(:),1,:),x(t(:),1,:),...
        x(tp1(:),1,:),x(tp2(:),1,:)],[],2);
    Nxmin=x./Nxmin;
    localmin(:,i,:) = ((Nxmin>=mcoef)+0);
    localmin(:,i,:) = localmin(:,i,:) .* (abs(x)>thresh);
    localmm=localmax+localmin;    
    localmm(localmm==2)=1;
end

end