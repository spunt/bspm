% make sure we are using the NIFTI toolbox



%addpath('/proj/braindata2/neurocinematics/toolboxes/NIFTI/')

%refimg='/proj/braindata2/neurocinematics/HarvardOxford/MNI152_T1_2mm_brain.nii';


function nii=fixOriginator(targetimg)

refimg='/Users/enrico/Documents/MATLAB/piazzolla/HarvardOxford/MNI152_T1_2mm_brain.nii';
%targetimg='clustered_rsa.nii';

ref=load_nii(refimg);

target=load_nii(targetimg);


target.hdr.dime.pixdim = ref.hdr.dime.pixdim;
target.hdr.dime.scl_slope = ref.hdr.dime.scl_slope;
target.hdr.dime.xyzt_units = ref.hdr.dime.xyzt_units;
target.hdr.hist = ref.hdr.hist;


target.original.hdr.dime.pixdim = ref.original.hdr.dime.pixdim;
target.original.hdr.dime.scl_slope = ref.original.hdr.dime.scl_slope;
target.original.hdr.dime.xyzt_units = ref.original.hdr.dime.xyzt_units;
target.original.hdr.hist = ref.original.hdr.hist;

%fileout = [target.fileprefix '_W_ORIGINATOR.nii'];

%save_nii(target,fileout);

nii=target;
