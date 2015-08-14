function vol = bramila_detrend(cfg)
% INPUT
%   cfg.infile = location where the subject NII file (4D)
%   cfg.vol = input can be a 4D volume, if this exists, then infile is not loaded but used as ID (optional)
%   cfg.write = 0 (default 0, set to 1 if you want to store the volume)
%   cfg.detrend_type = type of detrend, default is 'linear-demean', other options are 'spline' and 'linear-nodemean'
%   cfg.TR = TR (mandatory if detrend spline is used)
% OUTPUT
%   vol = a 4D volume detrended

if(isfield(cfg,'vol'))
	data=cfg.vol;
	% add check that it's a 4D vol
elseif(isfield(cfg,'infile'))
	nii=load_nii(cfg.infile);
	data=nii.img;
end
data=double(data);

type='linear-nodemean';
if(isfield(cfg,'detrend_type'))
    type=cfg.detrend_type;
end

% resize the data into a 2-dim matrix, time in first dimension
kk=size(data);
if(length(kk)==4)
    T=kk(4);
    tempdata=reshape(data,[],T);
    tempdata=tempdata';
    fprintf('Detrending data...');
else
    T=kk(1);
    tempdata=data;
end

m=mean(tempdata,1);
switch type
    case 'linear-demean'
        tempdata=detrend(tempdata);
    case 'linear-nodemean'
        tempdata=detrend(tempdata);
        for row=1:T
            tempdata(row,:)=tempdata(row,:)+m;
        end
    case 'spline'
        error('not implememented')
        % add here code
end

% resize the data back

if(length(kk)==4)
    tempdata=tempdata';
    vol=reshape(tempdata,kk);
    fprintf(' done\n');
else
    vol=tempdata;
end

if cfg.write==1 || nargout<1
    bramila_savevolume(cfg,vol,'EPI volume after detrending','mask_detrend.nii');
end
