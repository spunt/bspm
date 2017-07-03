This is a somewhat stand alone series of MATLAB scripts to compare modules (i.e. subnetworks or just set of voxels with a label) with known modules based on the literature.
The scripts are coded for images coregistered with FSL 2mm MNI 152 template. It can be easily adapted to other voxel sizes.

The current comparison is done against three reference atlases
1) Yeo et al. (2011) The organization of the human cerebral cortex estimated by intrinsic functional connectivity. J Neurophysiol 106(3):1125-65, 2011.
Data downloaded from http://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011
Data was coregistered to MNI_152_2mm template with FSL. Refer to FSL_yeo_MNI.sh script to replicate the results.

2a) Power et al. (2011) Functional network organization of the human brain. Neuron. 72(4):665-78. 
Data downloaded from http://www.nil.wustl.edu/labs/petersen/Resources.html File Consensus264.xls Columns G,H,I,AF,AK were stored in local Matlab file Power2011.mat
2b) Cole et al. (2013) Multi-task connectivity reveals flexible hubs for adaptive task control. Nature Neuroscience 16, 1348–1355 (2013)
Data obtained from supplementary materials of the paper http://www.nature.com/neuro/journal/v16/n9/extref/nn.3470-S1.pdf Table S9. This is the same as 2a with the difference that module with ID 2,6 and 13 were not taken into consideration
For both 2a) and 2b) volumes where created as spheres with 1cm diameter centered in the specified locations.

3) Gordon et al. (2014) Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations. Cerebral Cortex. 
Data downloaded from http://www.nil.wustl.edu/labs/petersen/Resources.html File Parcels.zip, file parcels.xls for the labels and file Parcels_MNI_222.nii for coordinates.

There are two scores returned per each reference atlas for each module in the input image
i. The Pearson's correlation between the spatial maps. Since number of voxels in each reference networks is different, only voxels with an actual label are considered.
ii. The ratio of overlap computed as the percentage of voxels in a certain given reference module

The Pearson's correlation seems to work quite well and has been used before in the literature (Smith 2009 PNAS).

