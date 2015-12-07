function [] = WaveletDespike(Inii,Opref,varargin)
%
% FUNCTION:     WaveletDespike -- Wavelet Despikes time series, as
%                                 described in:
%
% Patel, AX. et al (2014). A wavelet method for modeling and despiking 
% motion artifacts from resting-state fMRI time serires. NeuroImage.
%
%
% USAGE:        WaveletDespike(Inii,Opref,varargin)
%
% Inputs:       Inii      -- NIfTI file (3D+t dataset) containing time
%                            series to despike.
%               Opref     -- Output prefix for despiked NIfTI files and
%                            Spike Percentage (if specified).
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
%                            chains. Input must be a string containing one 
%                            of the following: 'conservative', 'moderate',
%                            'harsh'. [Default='moderate'].
%               nscale    -- Method for computing number of scales. Input
%                            must be a string containing one of the
%                            following: 'conservative','liberal','extreme'.
%                            [Default: 'liberal'].
%               compress  -- Binary flag indicating whether to compress out
%                            non-brain regions, 1, or not, 0, from input 
%                            NIfTI file before Wavelet Despiking. This 
%                            saves RAM and reduces runtime. [Default=1].
%               sp        -- Binary flag indicating whether to output the
%                            Spike Percentage for the dataset.
%                            [Default=1].
%               LimitRAM  -- Specify an upper bound of RAM usage (in Giga-
%                            bytes). Default is to use all available RAM.
%               verbose   -- Binary flag indicating whether to display
%                            incremental output from the algorithm, 1, or
%                            not, 0.
% 
% Outputs:      This function will write the following files to the current
%               directory:
%
%               Opref_wds.nii.gz   - Wavelet Despiked time series.
%               Opref_noise.nii.gz - Noise components removed in Wavelet
%                                    Despiking.
%               Opref_SP.txt       - The Spike Percentage (if specified at
%                                    input).
%
% EXAMPLE:      WaveletDespike('rest.nii.gz','rest','wavelet','la8',...
%                   'LimitRAM',5)
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

% ID: wavcorr.m 19 07-02-2014 BWTv1.1 axpatel


%% check nopts / parse extra opts

fname=mfilename;
if nargin<1
   help(fname); return;
end

if chkninput(nargin,[2,20],nargout,[0,0],fname)>=1
   return
end    

DefaultOpts=struct('wavelet','d4','threshold',10,'boundary','reflection',...
    'chsearch','moderate','nscale','liberal','compress',1,'SP',1,...
    'LimitRAM',0,'verbose',1);
Opts=parseInOpts(DefaultOpts,varargin);

%% check inputs

wavelets={'d4','d6','d8','d10','d12','d14','d16','d18','d20',...
      'la8','la10','la12','la14','la16','la18','la20',...
      'bl14','bl18','bl20','c6','c12','c18','c24','haar'};
boundaries={'periodic','circular','reflection'};
ok_chtype={'conservative','moderate','harsh'};
scaleopts={'conservative','liberal','extreme'};

err=struct();
err(1).inp=chkintype(Inii,'char',fname);
err(2).inp=chkintype(Opref,'char',fname);
err(3).inp=chkintype(Opts.wavelet,'char',fname,wavelets);
err(4).inp=chkintype(Opts.threshold,'numeric',fname);
err(5).inp=chkintype(Opts.boundary,'char',fname,boundaries);
err(6).inp=chkintype(Opts.chsearch,'char',fname,ok_chtype);
err(7).inp=chkintype(Opts.nscale,'char',fname,scaleopts);
err(8).inp=chkintype(Opts.compress,'numeric',fname,{'0','1'});
err(9).inp=chkintype(Opts.SP,'numeric',fname,{'0','1'});  
err(10).inp=chkintype(Opts.LimitRAM,'numeric',fname);
err(11).inp=chkintype(Opts.verbose,'numeric',fname',{'0','1'});
if sum(cat(1,err.inp))>=1;
   return
end

%% load NIfTI, parse MemorySolver for BatchMode

t=cputime;

[ts,Info,error]=ParseInNii(Inii,'compress',Opts.compress);

if error==1 
    return
end
header;

[NX,NTS]=size(ts);

if strcmpi(Opts.boundary,'reflection');
    NXref=NX*2;
    [BatchMode,freeRAM]=MemorySolver(NTS,NXref,'LimitRAM',Opts.LimitRAM);
else 
    [BatchMode,freeRAM]=MemorySolver(NTS,NX,'LimitRAM',Opts.LimitRAM);
end

if isempty(BatchMode);
    return; 
end

%% run wavelet despiking in normal mode

if BatchMode==0
    if Opts.SP==1
        [clean,noise,SP]=wdscore(ts,'wavelet',Opts.wavelet,'threshold',...
            Opts.threshold,'boundary',Opts.boundary,'chsearch',...
            Opts.chsearch,'nscale',Opts.nscale,'verbose',Opts.verbose);
    else
        [clean,noise]=wdscore(ts,'wavelet',Opts.wavelet,'threshold',...
            Opts.threshold,'boundary',Opts.boundary,'chsearch',...
            Opts.chsearch,'nscale',Opts.nscale,'verbose',Opts.verbose);
    end
end

%% run wavelet despiking in batch mode

if BatchMode>=1
    
    Ind=ParseInds(NTS,BatchMode);
    Batchstr=num2str(BatchMode);
    clean=zeros(NX,NTS);
    noise=zeros(NX,NTS);
    
    fprintf('Initialising Batch Mode ...\n')
	
    if Opts.LimitRAM==0
        fprintf('Using all available RAM. WARNING: your computer \n')
        fprintf('may be slow whilst this program is running.\n')
    else
        fprintf('Capping RAM usage at %s Gb ...\n',...
            num2str(freeRAM));
    end
        
    for ii=1:BatchMode
        
        fprintf('\nDespiking batch %s of %s.\n', num2str(ii),...
            Batchstr);
        
        tts=ts(:,Ind(1,ii):Ind(2,ii));
        
        if Opts.SP==1
            [clean(:,Ind(1,ii):Ind(2,ii)),noise(:,Ind(1,ii):Ind(2,ii)),...
                SP(:,ii)]=wdscore(tts,'wavelet',Opts.wavelet,...
                'threshold',Opts.threshold,'boundary',Opts.boundary,...
                'chsearch',Opts.chsearch,'nscale',Opts.nscale,'verbose',...
                Opts.verbose);
            SP(:,ii)=SP(:,ii).*((Ind(2,ii)-Ind(1,ii)+1)/NTS);
        else
            [clean(:,Ind(1,ii):Ind(2,ii)),noise(:,Ind(1,ii):Ind(2,ii))]...
                =wdscore(tts,'wavelet',Opts.wavelet,'threshold',...
                Opts.threshold,'boundary',Opts.boundary,'chsearch',...
                Opts.chsearch,'nscale',Opts.nscale,'verbose',Opts.verbose);
        end     
    end
    
    if exist('SP','var')
        SP=sum(SP,2);
    end
end
    
%% Write files

WriteOutNii(clean,strcat(Opref,'_wds.nii.gz'),Info);
WriteOutNii(noise,strcat(Opref,'_noise.nii.gz'),Info);

if exist('SP','var')
    dlmwrite(sprintf('%s_SP.txt',Opref),SP,' ');
end

footer(t);

end