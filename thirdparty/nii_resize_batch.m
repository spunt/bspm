function nii_resize_batch;
% Resize many images
%
%YOU CAN EDIT SUBJ TO BE LIST OF ALL IMAGES...
%  subj = strvcat('/usr/a/im1.nii','/usr/a/im2.nii','/usr/a/im3.nii');


%select files with a dialog
    subj = spm_select(inf,'image','Select images to compute resize');

 for i=1:size(subj,1)

   nii_resize_img(deblank(subj(i,:)),[2 2 2],[  -78 -112  -50;  78   76  85]);

end;
