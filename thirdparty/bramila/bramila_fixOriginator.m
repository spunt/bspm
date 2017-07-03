function nii=bramila_fixOriginator(targetimg)
% Usage:
%	nii = bramila_fixOriginator(filename);
%

% refimg could be improved
% it could have option to save already
refimg='/m/nbe/scratch/braindata/shared/toolboxes/HarvardOxford/MNI152_T1_2mm_brain.nii';

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
