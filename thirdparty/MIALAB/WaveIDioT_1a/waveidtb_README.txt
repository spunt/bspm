WaveIDioT toolbox Version 1a.

WaveIDioT is a Matlab toolbox allowing for improved 3-D denoising 
of fMRI data sets using a wavelet-based hierarchical approach. With 
the help of this toolbox, the user can easily smooth the data without
losing essential spatial details (edges, shape of gyri etc.) right 
after the other preprocessing steps have been applied. The GUI provides
an easy-to-use tool for beginner's whereas the batch scripts will soon 
be available for experienced users. 
Please send comments to vcalhoun@mrn.org or skhullar@mrn.org

Authors: Siddharth Khullar (skhullar@mrn.org)
		Chester F. Carlson Center for Imaging Science
		Rochester Institute of Technology, Rochester, NY, USA

	   Vince Calhoun, Ph.D (vcalhoun@mrn.org)
		Mind Research Network
		Albuquerque, NM, USA

Copyright 2011 Siddharth Khullar, Vince Calhoun.

References:

1. S. Khullar et al. , Wavelet-based fMRI analysis: 3-D denoising, signal separation and validation metrics. 
	Neuroimage, vol. 54, no. 4, Feb 2011.

2. S. Khullar et al. , Wavelet-based denoising and independent component analysis for improving multi-group
	inference in fMRI data. IEEE's International Symposium Biomedical Imaging, Chicago, IL, USA, March 2011.


Note: Requires MATLAB Version 7 or higher.
	SPM5/8 should be installed and added to the path in MATLAB
	in order to use this toolbox.

**************************************************************
Steps to use WAVEIDT toolbox for denoising fMRI data:

SETUP DENOISING FOR THE FIRST TIME.

1. After downloading the toolbox. Please unzip the directory and add the directory using "Set Path" or "addpath()" command.

Make sure Step 1 is completed (and SPM5/8 installed) otherwise it will not be possible to proceed forward.

2. Type WaveIDT1 on the command line of MATLAB (Case sensitive).

3. After the GUI opens, hit "SETUP DENOISING". 
   A new dialog box containing the directory structure will open.
   Please select a directory (by browsing) where the Denoising parameters MAT file 
   will be stored and accessed for later use. Press OK.

Note: If you do not select an output directory in Step 3, the present working directory 
	will be treated as an output directory.

4. A new window opens - "SPECIFY PARAMETERS". Here we specify the parameters required for denoising.
	
	- Type of Wavelet basis: Our work in [1] and [2] uses sym2 as the default as done here as well. 
					 Other options are available as a drop down menu and can be selected by 
					 the users if they wish.
					 Default: sym2 (as in MATLAB).

	- Number of decomposition levels : This is to set the number of wavelet decomposition levels used
							for hierarchical denoising explained in [1] and [2]. We suggest
							to set it to 2 or 3 for efficient computation, faster processing.
							Default: 2 levels

	- Select Data : Press Select, and chose the method of selection from the dropdown menu.
			    
			    SPM_SELECT_GUI: Opens the SPM5's generic window to select fMRI files.
			    It is suggested to use SPM's generic search strings (and recursive properties)
			    to find the files efficiently in case they are stored in a directory structure,
			    such as studies/subjects/sessions.

			    DATA_IN_ONE_FOLDER: If the user knows how the data is organized, this option is suggested.
							The following options will help user find the data recursively by 
							specifying the following:

						-- PREFIX of the data: represents the file pattern to be searched 
						                       for(has to one or more characters).

						-- Data FORMAT: What is format images files are in 
						   NIFTI(.nii) or ANALYZE (.img). Please use SPM_SELECT_GUI if you wish
       					   to select files with different formats.

						-- Select the root search directory: Select the root directory under which
						   all data is stored. This is "ONE" folder to be searched recursively for 
						   image files having the same pattern (PREFIX) and FORMAT(.nii or .img).

	After the Progress bar closes, please press DONE.				   

5. Press Run Denoising and follow the progress presented through the command window.

USE EXISTING PARAMETERS (saved as mat file) - The "Load Parameters" button.

1. Select the MAT file using the dialog box. The filename is fixed as:

waveidtb_ParameterInfo_<DateStored>_<TimeStored(HHMMSS)>.mat and stored in the output 
folder the user selected in the Step 1 before.

2. Press Run Denoising and follow the progress presented through the command window.


If you have any questions, please email skhullar@mrn.org.





