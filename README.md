bspm
====

This is the code I use to preprocess and analyze fMRI data and mostly relies on the routines that ship with [SPM12](http://www.fil.ion.ucl.ac.uk/spm/). Some of the functions prepended with `bspm_` are wrappers for direct calls to the SPM functions (e.g., `bspm_checkreg` is mostly lifted from `spm_check_registration`), but most are designed to make it easier to build the `matlabbatch` structures for various modules in SPM's batching system. I made this public because, well, there's no reason for me not to make it public. But please proceed with caution as I cannot guarantee everything is on the up-and-up.

If you encounter an error that you cannot solve yourself, let me know and I will do my best to help you solve it. If you encounter an error that you _can_ solve yourself, then I'd be grateful if you'd let me know. 

Three cheers for open source:

- One.
- Two.
- Three.