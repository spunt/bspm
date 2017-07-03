export 'FSLOUTPUTTYPE'='NIFTI' 
infile=$1;

mni1mm=/share/apps/fsl/5.0.9/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz
mni2mm=/share/apps/fsl/5.0.9/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz
ref=$mni2mm;

echo flirt -in $infile -ref $ref -out $infile"_toFSL.nii" -init spm2fsl.mat -applyxfm
flirt -in $infile -ref $ref -out $infile"_toFSL.nii" -init spm2fsl.mat -applyxfm
