function stats=bramila_oneSamplePermTest(cfg)
% usage
%	stats = bramila_oneSamplePermTest(cfg)
%
%	INPUTS
%		cfg.data (this has priority over infiles/mask)
%			One column per subject
%		OR
%
%   	cfg.infiles = cell vector with single volumes for each subject
%		cfg.mask = analysis mask (or the brain MNI mask)
%
%		cfg.mean = the expected mean for your data to satisfy the exchangeability condition (Nichols2002) or your chance level
%		cfg.niter = number of permutations
%		cfg.limits = [0 1] for accuracies, [-1 1] for correlations, [-inf +inf] for tvalues
%
%		stats.mean = the group mean
%		stats.maxTH = thresholds from the max statistics at 0.05 0.01 0.001 0.0001
%		stats.minTH = same but for the left tail
%
%	TODO
%		 input checking


% check the limits
EC=0;
if(any(isinf(cfg.limits)))
	disp('you are testing for t values or other values that can go up to infinite')
else
	if(mean(cfg.limits)~=cfg.mean)
		disp('your mean is not the actual mean of the limits, we will do a reshaping of intervals')
		EC=1;
	else
		disp('your mean is also the mean of the limits, no reshape needed')
	end
end

if(isfield(cfg,'data') && ~isempty(cfg.data))
	% use data
	data=cfg.data;
	NS=size(data,2);
	kk=[size(data,1) 1];
	inmask=1:kk(1);
else
	NS=length(cfg.infiles);

	mask=load_nii(cfg.mask);
	inmask=find(mask.img>0);

	for s=1:NS
		disp(['loading subject ' num2str(s) ' ' cfg.infiles{s}])
		temp=load_nii(cfg.infiles{s});
		data(:,s)=temp.img(:);
		kk=size(temp.img);
	end
end

data(find(isnan(data)))=0;

temp=mean(data,2);
stats.mean=	reshape(temp,kk);

% set seed
for iter=1:cfg.niter
	% generate random vectors of ones and minus ones
	temp=sign(randn(NS,1));
	if(EC==0)
		surromean=(data-cfg.mean)*temp/length(temp)+cfg.mean;
	else
		surromean=0;
		for t=1:length(temp)
			if(temp(t)==1)
				surromean=surromean+data(:,t);
			else
				temp2=temp(t)*(data(:,t)-cfg.mean); % flipped and removed mean
				if(0)
				negids=find(temp2<0);
				posids=find(temp2>0);
				temp2(negids)=abs(cfg.limits(1)-cfg.mean)*temp2(negids)/abs(cfg.limits(2)-cfg.mean);
				temp2(posids)=abs(cfg.limits(2)-cfg.mean)*temp2(posids)/abs(cfg.limits(1)-cfg.mean);
				end
				temp2=temp2+cfg.mean;
				surromean=surromean+temp2;
				%save debug
			end
		end
		surromean=surromean/length(temp);
	end
	Mperms(iter,1)=max(surromean(inmask));
	mperms(iter,1)=min(surromean(inmask));
end
stats.maxTH = prctile(Mperms,100-[5 1 0.1 0.01]);
stats.minTH = prctile(mperms,[5 1 0.1 0.01]);

