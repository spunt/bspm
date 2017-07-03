clear all
close all

subnii=load_nii('/scratch/braindata/shared/toolboxes/HarvardOxford/HarvardOxford-sub-prob-2mm.nii');

% white matter is in volumes 1 and 12
wm=double(subnii.img(:,:,:,1))/100+double(subnii.img(:,:,:,12))/100;
save_nii(make_nii(wm,[2 2 2]),'HO/white.nii')
temp(:,:,:,1)=wm;

% csf, or to be more precise "Ventricles", is in volumes 3 and 14
csf=double(subnii.img(:,:,:,3))/100+double(subnii.img(:,:,:,14))/100;
save_nii(make_nii(csf,[2 2 2]),'HO/csf.nii')
temp(:,:,:,2)=csf;

% cortical grey matter is in volumes 2 and 13
gr=double(subnii.img(:,:,:,2))/100+double(subnii.img(:,:,:,13))/100;
save_nii(make_nii(gr,[2 2 2]),'HO/grey_co.nii')
temp(:,:,:,3)=gr;

% subcortex
tempsub=subnii.img;
tempsub(:,:,:,[1 12 2 13 3 14])=[]; % we get rid of the stuff already processed, what is left is the subcortex probabilities
sub=sum(double(tempsub),4)/100;
save_nii(make_nii(sub,[2 2 2]),'HO/grey_sub.nii')
temp(:,:,:,4)=sub;

% we now do it for the cerebellum
cenii=load_nii('/scratch/braindata/shared/toolboxes/HarvardOxford/Cerebellum-MNIfnirt-prob-2mm.nii');
ce=sum(double(cenii.img),4)/100;
save_nii(make_nii(ce,[2 2 2]),'HO/grey_ce.nii')
temp(:,:,:,5)=ce;

% we now make it non-overlapping, by taking the maximum prob at each voxel
outmask=zeros(size(temp));
for x=1:size(temp,1)
	disp(num2str(x))
	for y=1:size(temp,2)
		for z=1:size(temp,3)
			thisval=squeeze(temp(x,y,z,:));
			if(sum(thisval)==0)
				continue;
			end
			M=max(thisval);
			if(length(M)>1)
				error('same max')
			end
			outmask(x,y,z)=find(M==thisval);
		end
	end
end
save_nii(make_nii(outmask,[2 2 2]),'HO/mask.nii')


