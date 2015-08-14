function Nimg_vec = waveidtb_neighborhood3d(img,nhood) %#codegen
%  NEIGHBORHOOD3D
% Is equivalent to colfilt operation in 3D
%
%   This function computes the neighborhood in 3D. NHOOD can
%   be 27-connected only.
%   The output image Nimg_vec will have the size of  r*c*d x nhood
%   where r, c and d are the number of rows, columns and depths


[R,C,Z]=size(img);
Nimg = ones(R,C,Z,nhood,'single');
r0=2:R-1;
c0=2:C-1;
z0=2:Z-1;

Nimg(r0,c0,z0,14) = img(r0,c0,z0);
Nimg(r0,c0,z0,10) = img(r0-1,c0-1,z0);
Nimg(r0,c0,z0,13) = img(r0-1,c0,z0);
Nimg(r0,c0,z0,16) = img(r0-1,c0+1,z0);
Nimg(r0,c0,z0,17) = img(r0,c0+1,z0);
Nimg(r0,c0,z0,18) = img(r0+1,c0+1,z0);
Nimg(r0,c0,z0,15) = img(r0+1,c0,z0);
Nimg(r0,c0,z0,12) = img(r0+1,c0-1,z0);
Nimg(r0,c0,z0,11) = img(r0,c0-1,z0);
% five voxels on top plane
Nimg(r0,c0,z0,5) = img(r0,c0,z0-1); % center
% 4-neighbors
Nimg(r0,c0,z0,4)= img(r0-1,c0,z0-1);
Nimg(r0,c0,z0,8)= img(r0,c0+1,z0-1);
Nimg(r0,c0,z0,6)= img(r0+1,c0,z0-1);
Nimg(r0,c0,z0,2)= img(r0,c0-1,z0-1);
% five voxels on bottom plane
Nimg(r0,c0,z0,23)= img(r0,c0,z0+1); % center
% 4-neighbors
Nimg(r0,c0,z0,22)= img(r0-1,c0,z0+1);
Nimg(r0,c0,z0,26)= img(r0,c0+1,z0+1);
Nimg(r0,c0,z0,24)= img(r0+1,c0,z0+1);
Nimg(r0,c0,z0,20)= img(r0,c0-1,z0+1);
% four corners-top
Nimg(r0,c0,z0,1) = img(r0-1,c0-1,z0-1);
Nimg(r0,c0,z0,7) = img(r0-1,c0+1,z0-1);
Nimg(r0,c0,z0,9) = img(r0+1,c0+1,z0-1);
Nimg(r0,c0,z0,3) = img(r0+1,c0-1,z0-1);
% four corners bottom
Nimg(r0,c0,z0,19) = img(r0-1,c0-1,z0+1);
Nimg(r0,c0,z0,25) = img(r0-1,c0+1,z0+1);
Nimg(r0,c0,z0,27) = img(r0+1,c0+1,z0+1);
Nimg(r0,c0,z0,21) = img(r0+1,c0-1,z0+1);

Nimg_vec = reshape(Nimg,R*C*Z, 27)';
