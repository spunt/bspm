function pow=bramila_ts_power(cfg)
	% takes a config struct with fields
	% cfg.niifile = path to a valid nifti file
	% cfg.toi=times of interest (interesting if you need to discard some time points)

	% need to implemant variable checks and throw errors
	nii=load_nii(cfg.niifile);
	img=double(nii.img);
	T=size(img,4);
	% now toi is forced to all time points
	toi=1:T;

	temp=reshape(img,[],T);
	temp=temp';	% time on the first dimension
	temp=temp(toi,:);
	pow=sum(temp.^2,1)/T;
	pow=reshape(pow,size(img,1),size(img,2),size(img,3));

