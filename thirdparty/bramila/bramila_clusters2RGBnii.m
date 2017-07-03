function bramila_clusters2RGBnii(cfg)
%
%	cfg.infile = path to nii file
%	cfg.colors = colormap file
% cfg.output_prefix= a string such './niis/modules' so that it will produce three files named ./niis/modules_R.nii ./niis/modules_G.nii ./niis/modules_B.nii
if(isfield(cfg,'output_prefix')==0)
	error('You need to specify output_prefix');
end

if(isa(cfg.output_prefix,'char')==0 || isempty(cfg.output_prefix) || strcmp(cfg.output_prefix,''))
	error('output_prefix cannot be empty and must be a string')
end



nii=load_nii(cfg.infile);
nii.img=round(nii.img);
ud=unique(nii.img);
%ud(find(ud)==0)=[];


if(any(ud<0)) error('only positive cluster IDs'); end

if(isfield(cfg,'colors')==0);
	colors=cbrewer('qual','Set1',length(ud));
else
	colors=cfg.colors;
end

colors=[0 0 0;colors]; % adding the black for the zero

% interpolating spaces between modules
allimg=zeros([size(nii.img) length(ud)]);
for uID=1:length(ud)
	tempvol=zeros(size(nii.img));
	tempvol(find(nii.img==ud(uID)))=1;
	smo=5;
	tempvol=smooth3(tempvol,'gaussian',[ smo smo smo]);
	allimg(:,:,:,uID)=tempvol;
end

% find the best match

out=zeros(91,109,91);
    for x=1:91
        disp(num2str(x))
        for y=1:109
            for z=1:91
                v=squeeze(allimg(x,y,z,:));

                if(sum(v)==0)
                    out(x,y,z)=0;
                else
                    v=v+eps*randn(size(v));
                    M=max(v);
                    temp=find(v==M(1)); %% add code that in case of 2 top vals, it finds what the non interpolated image would put...
                    out(x,y,z)=ud(temp(1));
                end
            end
        end
    end

save_nii(make_nii(out),'debug.nii')
% note that at the moment is going to give error if the cluster IDs are not increasing integers


RGBvol=zeros([91 109 91 3]);
    for uID=1:length(ud)
        ids=find(out==ud(uID));
        color=colors(uID,:);
        for rgb=1:3
            temp=RGBvol(:,:,:,rgb);
            temp(ids)=color(rgb);
            RGBvol(:,:,:,rgb)=temp;
        end
    end

    str={'R', 'G', 'B'};

    for rgb=1:3
        cfg.filename=[cfg.output_prefix '_' str{rgb} '.nii'];
		disp(['Writing ' cfg.filename]);
        save_nii(make_nii(RGBvol(:,:,:,rgb),[2 2 2]),cfg.filename);

        nii=bramila_fixOriginator(cfg.filename);
        save_nii(nii,cfg.filename);
    end

