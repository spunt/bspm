function matlabbatch = bspm_realign(epi_images, resliceflag)
% BSPM_REALIGN
%
% USAGE: bspm_realign(epi_images, resliceflag)
%
% If resliceflag = 1, epi images will be re-sliced
%

% ---------------- Copyright (C) 2014 ----------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1
   disp('Give me arguments! bspm_realign(epi_images)');
   return
elseif nargin<2
    resliceflag=0;
end

resliceflag = resliceflag*2;

if ischar(epi_images{1})
    nsess = 1;
    charflag = 1;
else
    nsess = length(epi_images);
    charflag = 0;
end

for s = 1:nsess

    if charflag
        c_epi_images = epi_images;
    else
        c_epi_images = epi_images{s};
    end

    % fix end of image filename cell array
    for i = 1:length(c_epi_images);
        c_epi_images(i) = cellstr([c_epi_images{i}, ',1']);
    end

    % put images into job
    matlabbatch{1}.spm.spatial.realign.estwrite.data{s} = c_epi_images;

end

matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {};
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [resliceflag 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = '';

% run job
if nargout==0,  spm_jobman('initcfg'); spm_jobman('run',matlabbatch); end


end
