bspm
====

This is the code I use to preprocess and analyze fMRI data and mostly relies on the routines that ship with [SPM12](http://www.fil.ion.ucl.ac.uk/spm/). Some of the functions prepended with `bspm_` are wrappers for direct calls to the SPM functions (e.g., `bspm_checkreg` is mostly lifted from `spm_check_registration`), but most are designed to make it easier to build the `matlabbatch` structures for various modules in SPM's batching system. 

I figured I'd make it open source because, well, I believe in open source and, in fact, much of what is in this repository is taken from other people's open source software. But because this hasn't been written _for_ the public, I can't gaurantee everything will work on your machine. But, I _can_ gaurantee a higher likelihood of working if you make sure that all folders and subfolders other than `thirdparty` is in your MATLAB path when using it. As for the `thirdparty` software, for most functions you'll need to add a version of `spm` to your path if you don't already have your own install. If you run into an undefined function error with any of the `bspm` functions, then you may be able to quickly solve it by finding the relevant function in the `thirdparty` folder and adding the toolbox it is in to your path. 

If you encounter an error that you cannot solve yourself, let me know and I will do my best to help you solve it. If you encounter an error that you _can_ solve yourself, then I'd be grateful if you'd let me know. 

Three cheers for open source:

- One.
- Two.
- Three.