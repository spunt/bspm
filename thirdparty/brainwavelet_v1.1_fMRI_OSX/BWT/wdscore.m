function [clean,noise,SP] = wdscore(ts,varargin)
%
% FUNCTION:     wdscore  -- Wavelet Despiking core algorithm, as described
%                           in:
%
% Patel, AX. et al (2014). A wavelet method for modeling and despiking 
% motion artifacts from resting-state fMRI time serires. NeuroImage.
%
%
% USAGE:        wdscore(ts,varargin)
%
% Inputs:       ts        -- Input matrix of time series with dimensions 
%                            NX x Nts, where NX is the number of time 
%                            points, and Nts is the number of time series
%                            to despike.
%
%               Additional Input Options:
%               (These must be specific as MATLAB string-value pairs.)
%
%               wavelet   -- Wavelet to use for wavelet transform. Input
%                            must be a string containing one of the 
%                            following:
%                            'd4','d6','d8','d10','d12','d14','d16','d18',
%                            'd20','la8','la10','la12','la14','la16',
%                            'la18','la20','bl14','bl18','bl20','c6','c12',
%                            'c18','c24','haar'. [Default='d4'].
%               threshold -- Threshold for maximal and minimal wavelet
%                            coefficients. [Default=10].
%               boundary  -- Boundary condition to apply. Input must be a
%                            string containing one of the following: 
%                            'circular','reflection'.
%                            [Default='reflection'].
%               chsearch  -- Rules for identifying maxima and minima 
%                            chains. [Default='moderate'].
%               nscale    -- Method for computing number of scales. Input
%                            must be a string containing one of the
%                            following: 'conservative','liberal','extreme'.
%                            [Default: 'liberal'].
%               verbose   -- Binary flag indicating whether to display
%                            incremental output from the algorithm [1] or
%                            not [0].
% 
% Outputs:      clean     -- Matrix of Wavelet Despiked time series, with
%                            dimensions NX x Nts.
%               noise     -- Matrix of recomposed signals removed by the
%                            Wavelet Despiking algorithm, with dimensions
%                            NX x Nts.
%               SP        -- The Spike Percentage for each time point.
%
%
% EXAMPLE:      [clean,noise,sp]=wdscore(ts,'wavelet','la8')
%
% AUTHOR:       Ameera X Patel
% CREATED:      26-12-2013
%   
% CREATED IN:   MATLAB 7.13
%
% REVISION:     19
%
% COPYRIGHT:    Ameera X Patel (2013, 2014), University of Cambridge
%
% TOOLBOX:      BrainWavelet Toolbox v1.1

% ID: wavcorr.m 19 05-02-2014 BWTv1.1 axpatel


%% check nopts / parse extra opts

fname=mfilename;
if nargin<1
    help(fname); return;
end

if chkninput(nargin,[1,13],nargout,[2,3],fname)>=1
    clean=[]; noise=[];
    if nargout==3; SP=[]; end
    return
end

DefaultOpts=struct('wavelet','d4','threshold',10,'boundary','periodic',...
    'chsearch','moderate','nscale','liberal','verbose',1);
Opts=parseInOpts(DefaultOpts,varargin);

if chkintype(Opts.verbose,'numeric',fname,{'0','1'})>=1;
    clean=[]; noise=[];
    if nargout==3; SP=[]; end
    return;
end

if Opts.verbose==1
    ProgMsg('inp',nargout)
end

%% check inputs

wavelets={'d4','d6','d8','d10','d12','d14','d16','d18','d20',...
      'la8','la10','la12','la14','la16','la18','la20',...
      'bl14','bl18','bl20','c6','c12','c18','c24','haar'};
boundaries={'periodic','circular','reflection'};
ok_chtype={'conservative','moderate','harsh'};
scaleopts={'conservative','liberal','extreme'};

err=struct();
err(1).inp=chkintype(ts,'double',fname);
err(2).inp=chkintype(Opts.wavelet,'char',fname,wavelets);
err(3).inp=chkintype(Opts.threshold,'numeric',fname);
err(4).inp=chkintype(Opts.boundary,'char',fname,boundaries);
err(5).inp=chkintype(Opts.chsearch,'char',fname,ok_chtype);
err(6).inp=chkintype(Opts.nscale,'char',fname,scaleopts);

if sum(cat(1,err.inp))>=1;
    clean=[]; noise=[]; 
    if nargout==3; SP=[]; end
    return
end

if Opts.verbose==1
    ProgMsg('done',nargout)
end

%% compute scales, initiate variables, and do 3-D wavelet transform

if Opts.verbose==1
    ProgMsg('modwt',nargout)
end

[NX,Nts]=size(ts);

if NX==1 && Nts==1;
    cprintf('_err','*** BrainWavelet Error ***\n')
    cprintf('err','Your time series is 1 time point long!\n',fname);
    clean=[]; noise=[];
    if nargout==3; SP=[]; end
    return
elseif NX==1 && Nts>1;
    ts=ts(:);
    NXc=NX; NX=Nts; Nts=NXc; clear NXc;
end

NJ=modwt_scales(NX,Opts.nscale,Opts.wavelet);
[m,sc,inf]=modwt(ts,Opts.wavelet,NJ,Opts.boundary,'RetainVJ',1);

if Opts.verbose==1
    ProgMsg('done',nargout)
end

%% check Nscales and chain search opts against scales

if NJ==1
    cprintf('err','\nCannot Wavelet Despike with only 1\n');
    cprintf('err','wavelet scale. Exiting...\n');
    clean=[]; noise=[];
    if nargout==3; SP=[]; end
    return
end

if (strcmpi(Opts.chsearch,'conservative') || ...
        strcmpi(Opts.chsearch,'moderate'))...
        && NJ==2
    cprintf('\nNB: %s chain search method with %s scales\n',...
        Opts.chsearch,num2str(NJ))
    cprintf('is equivalent to use of the harsh method');
    cprintf('Using harsh method for speed\n\n');
    Opts.chsearch='harsh';
end

%% compute phase-delay properties and apply cshift

if Opts.verbose==1
    ProgMsg('talc',nargout)
end
cshift=calccshift(Opts.wavelet,NJ);

NX=size(m,1);

% amend any cshifts greater than lts - all coefficients will be boundary.
iter=abs(min(ceil(cshift/NX)));
for i=1:iter
    ind=find(cshift<-NX);
    cshift(ind)=cshift(ind)+NX;
end
csmin=abs(cshift);
mcs=zeros(NX,NJ,Nts);

% circularly shift 3-D matrix of 1/2/3 x 4 x 5-D - use logical indexing.
for i=1:NJ
    mcs(:,i,:)=[m((1+csmin(i)):NX,i,:); m(1:csmin(i),i,:)];
end

if Opts.verbose==1
    ProgMsg('done',nargout)
end

%% find local maxima and minima

if Opts.verbose==1
    ProgMsg('lmm',nargout)
end

localmm=maxcoeff(mcs,Opts.threshold,0.5);
pmm=(localmm.*mcs>0)+0;
nmm=(localmm.*mcs<0)+0;

if Opts.verbose==1
    ProgMsg('done',nargout)
end

wl=5;  wlh=(wl-1)/2;
wh=3;  whh=(wh-1)/2;

%% parse denoising opts and find neighbourhood matrices

pmmp=circularize(pmm,'y',wlh);
nmmp=circularize(nmm,'y',wlh);
pmmp=zpad(pmmp,'x',whh);
nmmp=zpad(nmmp,'x',whh);
pnghbr=zeros(NX,NJ,Nts);
nnghbr=zeros(NX,NJ,Nts);

if Opts.verbose==1
    ProgMsg('nghbr',nargout)
end

if strcmpi(Opts.chsearch,'harsh')
    pnghbr=nghbrcount(pmmp,2,1,NX,NJ,Nts,1);
    nnghbr=nghbrcount(nmmp,2,1,NX,NJ,Nts,1);    
elseif strcmpi(Opts.chsearch,'moderate')
    pnghbr=mmcs_mod(pmm,pmmp,2,1,1);
    nnghbr=mmcs_mod(nmm,nmmp,2,1,1);
elseif strcmpi(Opts.chsearch,'conservative');
    pnghbr=mmcs_cons(pmm,pmmp,2,1,1);
    nnghbr=mmcs_cons(nmm,nmmp,2,1,1);
end

if Opts.verbose==1
    ProgMsg('done',nargout)
end

%% dump workspace variables (memory)

clear wl wlh whh wh pmmp nmmp localmm  mcs

%% chain search

if Opts.verbose==1
    ProgMsg('ch',nargout)
end

pres=zeros(NX,NJ,Nts);
nres=zeros(NX,NJ,Nts);
pres((pnghbr+pmm)==2)=1;
nres((nnghbr+nmm)==2)=1;
mmc=pres+nres;

if Opts.verbose==1
    ProgMsg('done',nargout)
end

%% dump workspace variables (memory)

clear pnghbr nnghbr pmm nmm pres nres

%% apply inverse phase-delay cshift

if Opts.verbose==1
    ProgMsg('pdsh',nargout)
end

mmcs=zeros(NX,NJ,Nts);
for i=1:NJ
    mmcs(:,i,:)=[mmc(NX-csmin(i)+1:NX,i,:); mmc(1:NX-csmin(i),i,:)];
end

if Opts.verbose==1
    ProgMsg('done',nargout)
end

%% find noise and signal coefficients

if Opts.verbose==1
    ProgMsg('sn',nargout)
end

wnoise=mmcs.*m;
mmcsi=abs(mmcs-1);
wsignal=mmcsi.*m;

if Opts.verbose==1
    ProgMsg('done',nargout)
end

%% compute imodwt 

if Opts.verbose==1
    ProgMsg('imodwt',nargout)
end

clean=imodwt(wsignal,sc,inf);
noise=imodwt(wnoise,sc,inf);

if Opts.verbose==1
    ProgMsg('done',nargout)
end

%% compute Spike Percentage

if nargout>=3
    if Opts.verbose==1
        ProgMsg('sp',nargout)
    end
    
    SP=SpikePercentage(mmc,m);
    
    if strcmpi(Opts.boundary,'reflection')
        NXh=NX/2;
        SP=(SP(1:NXh)+flipdim(SP(NXh+1:NX),1))/2;
    end
    
    if Opts.verbose==1
        ProgMsg('done',nargout)
    end
end

end