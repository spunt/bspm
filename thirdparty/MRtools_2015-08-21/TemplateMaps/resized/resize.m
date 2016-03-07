h1 = spm_vol('/Users/aaron/MATLAB_Scripts/spm12/canonical/single_subj_T1.nii');
h1.mat
h1.dt = [16 0];
h1.pinfo = [1 0 352]';
cd /Users/aaron/Dropbox/AaronsScripts/MRtools/TemplateMaps/
mkdir resized

[a b] = dir_wfp('*.nii');

for ii = 1:numel(b);
    h2 = spm_vol(b{ii});
    h1.fname = ['resized/' h2.fname];
    m = resizeVol2(h2,h1,3);
    m(isnan(m))=0;
    
    %m = m./(range(m(:))/2);
    spm_write_vol(h1,m);
end

cd resized