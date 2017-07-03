function groupmask=bramila_groupmask(cfg)
	% cfg.indata contains a list of files for the whole group
	Nsubj=length(cfg.indata);
	cfg.template='/triton/becs/scratch/braindata/shared/HarvardOxford/MNI152_T1_2mm_brain.nii';
	templ=load_nii(cfg.template);
	inmask=find(templ.img>0);
	groupmask=ones(size(templ.img));
	for sub=1:Nsubj
		Pow=bramila_tspower(cfg.indata{sub});
		p=prctile(Pow(inmask),2);
        smask=zeros(size(Pow));
        smask(find(Pow>=p))=1;
        groupmask=groupmask.*smask;
	end
