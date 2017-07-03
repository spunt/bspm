This folder contains tissue priors for the MNI152 2mm template size as per FSL templates, voxel dimensions of volume is 91 109 91.
The three volumes
	- grey.nii 	(grey matter)
	- csf.nii	(cerebral spinal fluid)
	- white.nii	(white matter)
are tissue priors (with probability values from 0 to 1) used for tissue regression. The tissue come from SPM8 'apriori' subfolder now moved to spm12b 'toolbox/fieldmap' subfolder.
There are other tissue priors in the HO (Harvard Oxford) subfolder, based on the harvard oxford probability atlas released with Fsl5.0. The matlab script HO_masks.m is the one used to generate these priors from the HO atlas files.


