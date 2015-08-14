function nii_fliplr (V);
% Generates left-right mirrorred image
%   nii_mirror('C:\dir\img.nii');
% You can also pass multiple images - each will be mirrored
% Example - scans from participant 1 and 2
%  nii_mirror(strvcat('C:\dir\p1.nii','C:\dir\p2.nii'));

if nargin <1 %no files specified
 V = spm_select(inf,'image','Select image to left-right flip');
end

disp('nii_fliplr: Flipping order of rows - please make sure this is the left-right axis');
for i=1:size(V,1)
 ref = deblank(V(i,:));
 [pth,nam,ext] = spm_fileparts(ref); 
 Vi  = spm_vol([ref,',1']);
 VO       = Vi; %create output header
 VO.fname = fullfile(pth,['RL' nam '.nii']);
 VO       = spm_create_vol(VO);
 xvect=zeros(1,Vi.dim(1));
 for i=1:Vi.dim(3), %flip every slice
    img      = spm_slice_vol(Vi,spm_matrix([0 0 i]),Vi.dim(1:2),0);
    for Yi=1:Vi.dim(2),
        px = (Yi-1)*Vi.dim(1);
        for Xi=1:Vi.dim(1),
            xvect(Xi) = img(px+Vi.dim(1)-Xi+1);
        end; %for X each column
        for Xi=1:Vi.dim(1),
            img(px+Xi) = xvect(Xi);
        end; %for X: each column

    end; %for Y: each row
    VO       = spm_write_plane(VO,img,i);
 end; %for Z %each slice
 
end; %for each image
return
