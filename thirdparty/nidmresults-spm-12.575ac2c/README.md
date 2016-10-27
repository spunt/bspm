
Test for NIDM-Results export in SPM
================

Testing procedures for SPM NIDM-Results export.

Test data is available at https://github.com/incf-nidash/nidmresults-examples/, you will need a local copy of this repository stored with [git lfs](https://git-lfs.github.com/):
```
git clone https://github.com/incf-nidash/nidmresults-examples.git
git lfs install
```

To run the test battery, follow those steps:
 1. Run the NIDM-Results export on your machine (within Matlab)
```
nidm_export_all('LOCAL_PATH_TO_NIDMRES_EX/nidmresults-examples/', 'LOCAL_PATH_TO_NIDMRES_SPM/nidm-results_spm/spmexport')
``` 
 2. Push the updates ttl and provn file to GitHub:
```
    git add -u
    git commit -m "description of updated feature"
    git push origin master
```
