fslswapdim Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz RL PA IS parc1mm.nii.gz
fslswapdim FSL_MNI152_FreeSurferConformed_1mm.nii.gz RL PA IS T1.nii.gz
flirt -in T1.nii.gz  -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii -out T12mm.nii.gz -omat one2two.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear
flirt -in parc1mm.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii -out YEO_parc2mm_fsl.nii.gz -applyxfm -init one2two.mat -interp trilinear
gunzip YEO_parc2mm_fsl.nii.gz
