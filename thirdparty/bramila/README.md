# BRAin and MInd LAb (BRAMILA) code repository

This is the code repository for the Brain and Mind Lab at the Department of Neuroscience and Biomedical Engineering at Aalto University. The main project stored in this repository is the BRAMILA fMRI preprocessing pipeline.

# BRAMILA fMRI preprocessing pipeline

## What is it
Enrico to add here

## How to use it
Link to docs on Aalto wiki

## How to report it in your paper
*This is party taken from a recently submitted paper, please avoid plagiarizing by rewriting and rephrasing*

Standard fMRI preprocessing steps were applied using the FSL software (www.fmrib.ox.ac.uk, version 5.0.9) and custom MATLAB code (BRAMILA pipeline v2.0, available at https://git.becs.aalto.fi/bml/bramila/). Briefly, EPI slices were firstly corrected for slice timing differences. Volumes were then corrected for head motion using MCFLIRT and coregistered to the Montreal Neurological Institute 152 2mm template in a two-step registration procedure using FLIRT: from EPI to participant’s anatomical image after brain extraction (9 degrees of freedom) and from anatomical to standard template (12 degrees of freedom). Further, spatial smoothing was applied with a Gaussian kernel of 6 mm full width at half maximum. To remove scanner drift, a 240-sec-long Savitzky-Golay filter ([Çukur, et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3929490/)) was applied. Furthermore, time series were also filtered with a high pass temporal filter at a cut-off frequency of 0.01 Hz. To further control for motion and physiological artefacts, BOLD time series were cleaned using 24 motion-related regressors, signal from deep white matter, ventricles and cerebral spinal fluid locations as described in [Power et al 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3849338/).



.
