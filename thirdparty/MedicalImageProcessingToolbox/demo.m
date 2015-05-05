%% Medical Image Processing Toolbox
% (c) Alberto Gomez 2013, King's College London
% alberto.gomez@kcl.ac.uk
% This software is provided as is. Please check the license file for
% details.

%% Basic image creation and visualization
load mri
orig = [0 0 0]';
sp = [1 1 1]';
orient = eye(3);
mri_im = ImageType(size(squeeze(D))', orig, sp, orient);
mri_im.data = squeeze(D);

figure,
mri_im.show()
axis equal;
view(40,30)
colormap(gray)

%% Extract an oblique slice

bounds = mri_im.GetBounds();
centroid =  (bounds([1 3 5])+bounds([2 4 6]))/2;
normal = [0.5 0.5 0.5]';
normal=normal/norm(normal);

slice = resliceImage(mri_im,'plane',normal,centroid,'interpolation','linear');

figure
mri_im.show();
hold on;
slice.show();
hold off;
colormap(gray)
axis equal;
view(40,30);

%% Downsample/resample the image

new_resolution = [2 2 2]';

mri_im_downsampled = resampleImage(mri_im,[],'spacing',new_resolution,'interpolation','linear');
figure,
subplot(1,2,1)
mri_im.show()
axis equal;
view(40,30)
title('Original')
subplot(1,2,2)
mri_im_downsampled.show()
axis equal;
view(40,30)
colormap(gray)
title('Downsampled')

