function cfg = bramila_filter(cfg)
% BRAMILA_FILTER - Bandpass filter fMRI time series
%	- Usage:
%	[vol, cfg] = bramila_filter(cfg) Returns the filtered volume and
%		filter specifications added to the cfg variable
%	- Input:
%	cfg is a struct with following parameters:
%   	cfg.infile = location where the subject NII file (4D)
%   	cfg.vol = input can be a 4D volume, if this exists, then infile is not loaded but used as ID (optional)
%		cfg.mask = mask volume (mandatory)
%   	cfg.write = 0 (default 0, set to 1 if you want to store the volume)
%   	cfg.TR = mandatory field in seconds
%		cfg.HPF = value in Hz, optional, default 0.01 Hz
%		cfg.LPF = value in Hz, optional, default 0.08 Hz
%		cfg.filtertype = 'butter' or 'fir' (default butter)
%		cfg.filterorder = only for butterworth filter, default 2
%	- Note:
%   We used the same filtering strategy as in Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048
if isfield(cfg,'vol') && ~isempty(cfg.vol)
	data=cfg.vol;
	% add check that it's a 4D vol
elseif isfield(cfg,'infile') && ~isempty(cfg.infile)
	nii=load_nii(cfg.infile);
	data=nii.img;
else
    error('No volume or infile found!') 
end


data=double(data);
mdata=mean(data,4); % mean in time;

if(~isfield(cfg,'TR'))
	error('cfg.TR is a mandatory field')
end
TR=cfg.TR;
Fs=1/TR;

% cut offs in Hz, default values [double check Power code]

A = [0 1 0];
DEV=[0.05 0.01 0.05];

F = [0,0.01,0.08,0.09];
if(isfield(cfg,'filter_limits'))
	%HPF=cfg.HPF;
    F = cfg.filter_limits;
end

HPF = F(2);
LPF = F(3);

if(isfield(cfg,'HPF'))
	HPF=cfg.HPF;
end
if(isfield(cfg,'LPF'))
	LPF=cfg.LPF;
end
if(LPF<=HPF)
	error('The freq for the low pass filter is smaller than the freq for the high pass filter')
end

F=[0,HPF, LPF,LPF+0.01]; % this should control that we don't go over nyquist

WRITE=0;
if(isfield(cfg,'write'))
	WRITE=cfg.write;
end

% FIR or BUTTER filters
FILTERTYPE='butter';
if(isfield(cfg,'filtertype'))
	FILTERTYPE=cfg.filtertype;
	% add check that it can only be 'butter' or 'fir' case insensitive
end

filterwithmatlab=1;
switch FILTERTYPE
	case 'butter'
		FILTERORDER=2;
		if(isfield(cfg,'filterorder'))
			FILTERORDER=cfg.filterorder;
		end
		
		hipasscutoff=HPF/(0.5/TR);
		lowpasscutoff=LPF/(0.5/TR);
		[b a]=butter(FILTERORDER,[hipasscutoff lowpasscutoff]);
		cfg.filter.butterfreq=[hipasscutoff lowpasscutoff];
		cfg.filter.butterorder=FILTERORDER;
	case 'fir'
        
		[N,Fo,Ao,W] = firpmord(F,A,DEV,Fs);
		if(mod(N,2)==1)
    		N=N+1;
		end
		% Design filter
		b=firpm(N,Fo,Ao,W);
		a=1;
		% store here in the cfg.filter other FIR specific parameters
		
		cfg.filter.N=N;
		cfg.filter.Fo=Fo;
		cfg.filter.Ao=Ao;
		cfg.filter.W=W;
	case 'fslhp'
		disp('Using FSL for high-pass temporal filtering');
		delete_EPI_temp=0;
		if(~isfield(cfg,'infile') || isempty(cfg.infile))
			fprintf('..creating tempfile for temporal filtering\n');
			cfg.infile = bramila_savevolume(cfg,cfg.vol,'temprorary file for temporal filtering','EPI_tempTfile.nii');
			delete_EPI_temp = 1;
		end
		infile=cfg.infile;
		outfile=[cfg.infile(1:end-4) '_HPF.nii' ];
		filterwithmatlab=0;
		sigm=round((1/HPF)/TR);
		setenv('FSLOUTPUTTYPE','NIFTI') % set output type to unarchived .nii
		command=['fslmaths ' infile ' -bptf ' num2str(sigm)  ' -1 ' outfile];
		disp(command)
		system(command);
		if delete_EPI_temp==1
			delete(cfg.infile);
            disp(['leaving  ' cfg.infile])
			rmfield(cfg,'infile');
		end
    otherwise
        error('Unknown filter type (only ''butter'', ''fir'' and ''fslhp''  allowed)');
end
if(filterwithmatlab)
    disp('Filtering with Matlab')
	cfg.filter.b=b;
	cfg.filter.a=a;
    
	% prepare the 4D data to be filtered
	siz=size(data);
	temp=reshape(data,[],siz(4));
	tsdata=double(temp');

	T = size(tsdata,1);
	m=mean(tsdata,1);	% storing the mean
	%m=repmat(m,T,1);
	%tsdata=tsdata-m;	% removing the mean
	for row = 1:T
		tsdata(row,:)=tsdata(row,:)-m;
	end

	mask=cfg.mask;	% this should be before since mask is mandatory
	maskID=find(mask>0);

	tsdataout=zeros(size(tsdata));

	fprintf('Filtering data...')
    if(strcmp(FILTERTYPE,'fir') && T<cfg.filter.N*3 && a==1 ) % FIR only
        disp('not using filtfilt')
        for v=1:size(tsdataout,2)
           if(ismember(v,maskID))
               temp=tsdata(:,v);
               frstpass=conv(b,temp);
               tsout=flipud(conv(b,flipud(frstpass)));
               tsout = tsout(length(b):end);
               tsout(T+1:end)=[];
               tsdataout(:,v)=tsout;
           end
        end
            
    else
        
        tsdataout(:,maskID)=filtfilt(b,a,tsdata(:,maskID));
    end
	fprintf(' done\n');

	% NOTE: Power code applies filter in bordered data (removing stuff at beginning and end): border data is better than no data, we keep it now but later for the FC measure we use only the inside data 
	%tsdataout=tsdataout+m; % add the mean back, useful for dvars computations
	for row = 1:T
		tsdataout(row,:)=tsdataout(row,:)+m;
	end
	tsdataout=tsdataout';
	vol=reshape(tsdataout,siz);
	cfg.vol = vol;

	if WRITE==1 || nargout<1
		cfg.infile = bramila_savevolume(cfg,vol,'masked, detrended, regressed and filtered EPI volume','mask_detrend_fullreg_filtered.nii');
		cfg.outfile = cfg.infile;
	end
else
	% if you are here, it's because you are using FSL to HP filter
	temp=load_nii(outfile);
	% we remove the mean and we readd th mean
	mtemp=mean(temp.img,4);
	for t=1:size(temp.img,4)
		temp.img(:,:,:,t)=temp.img(:,:,:,t)+mdata-mtemp;
	end
	cfg.vol=temp.img;
	cfg.infile = bramila_savevolume(cfg,cfg.vol,'masked, detrended, regressed and filtered EPI volume','mask_detrend_fullreg_filtered.nii');
	cfg.outfile=cfg.infile;
end





