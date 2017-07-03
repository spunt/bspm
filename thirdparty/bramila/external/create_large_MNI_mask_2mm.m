clear all;

nii_orig = load_nii('MNI152_T1_2mm_brain_mask.nii');
mask = nii_orig.img;
N = nnz(mask);

nii = load_nii('grey.nii');
grey = nii.img;
nii = load_nii('white.nii');
white = nii.img;
nii = load_nii('csf.nii');
csf = nii.img;

mask(grey>0.1)=1;
mask(csf>0.1)=1;
mask(white>0.1)=1;

NN=nnz(mask);

fprintf('old mask size %i, new mask size %i\n',N,NN);

nii_orig.img=mask;
save_nii(nii_orig,'MNI152_T1_2mm_brain_mask_ENLARGED.nii');

