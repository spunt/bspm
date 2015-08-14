GroupICATv4.0a Updates (June 21, 2015):

Image viewer and HTML report utilities are fixed to work on MATLAB R2008a. The following files are updated:
	a. icatb/gift.m
	b. icatb/icatb_helper_functions/icatb_image_viewer.m

GroupICATv4.0a Updates (June 17, 2015):

Option is provided to do paired t-test in stats on beta weights utility. File icatb/icatb_helper_functions/icatb_statistical_testing_TC.m is updated.


GroupICATv4.0a Updates (June 11, 2015):

1. MPOWIT (un-stacked) code is updated to use random initialization if STP initialization fails due to out of memory error. The following files are updated:
	a. icatb/icatb_analysis_functions/icatb_calculate_pca.p
	b. icatb/icatb_parallel_files/icatb_parCalculatePCA.p
2. Message is not printed to command window if verbose option is turned off. The following files are updated:
	a. icatb/icatb_analysis_functions/icatb_eig_symm.m
	b. icatb/icatb_analysis_functions/icatb_calculate_em_pca.m

GroupICATv4.0a Updates (May 05, 2015):

1. HTML viewer is fixed to handle components written in compressed format and/or written in analysis sub-directories. File icatb/icatb_helper_functions/icatb_gica_html_report.m is changed.
2. Error message "Cannot convert double value NaN to a handle" when using ICASSO on R2014b is fixed. File icatb/toolbox/icasso122/clusterhull.m is modified.

GroupICATv4.0a Updates (May 04, 2015):

1. Orthogonal viewer utility is not functional when accessed from display tools drop down box. Added case match for "orthogonal viewer" in 
icatb/icatb_helper_functions/icatb_utilities.m to fix the problem.

2. ICA parameter file is saved before executing MDL estimation tool. The following files are modifed:
	a. icatb/icatb_batch_files/icatb_read_batch_file.m
	b. icatb/icatb_helper_functions/icatb_estimateCompCallback.m
