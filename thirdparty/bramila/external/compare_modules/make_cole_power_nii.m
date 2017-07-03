clear all
close all
% Cole data was downloaded from
% http://www.colelab.org/cole-etal-2013/fullsupplementaryresults.pdf
%


 addpath('/Users/enrico/code/bramila/git/bramila/') % load bramila toolbox for MNI conversions / run once

files={
	'power2011_modules.csv'
	'cole2013_modules.csv'
	};

for f=1:length(files)
	data=load(files{f});

	outvol=zeros(91,109,91);
	for r=1:size(data,1)
		[x,y,z]=bramila_MNI2space(data(r,2),data(r,3),data(r,4));
		outvol(round(x),round(y),round(z))=data(r,1); % make sure x y z are integer indexes
		% let's put a sphere around 
		outvol(round(x)+[-1:1],round(y)+[-1:1],round(z)+[-1:1])=data(r,1);
		outvol(round(x)+[-2:2],round(y),round(z))=data(r,1);
		outvol(round(x),round(y)+[-2:2],round(z))=data(r,1);
		outvol(round(x),round(y),round(z)+[-2:2])=data(r,1);
    end
	filename=files{f};
	filename(end-3:end)=[];
	filename=[filename '.nii'];
	save_nii(make_nii(outvol),filename);
	nii=fixOriginator(filename);
	save_nii(nii,filename);
end

