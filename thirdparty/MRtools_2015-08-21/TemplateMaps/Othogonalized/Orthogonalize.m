cd /Users/aaron/Dropbox/AaronsScripts/MRtools/TemplateMaps
[pths files] = dir_wfp('/Users/aaron/Dropbox/AaronsScripts/MRtools/TemplateMaps/*.nii');

h = spm_vol(char(files));
m = FastRead(files);
% mm = FastRead(dir_wfp('/Users/aaron/Dropbox/AaronsScripts/MRtools/TemplateMaps/Othogonalized/*.nii'));
%%
mkdir Othogonalized/
delete Othogonalized/*.nii

vol = zeros(h(1).dim);
for ii = 1:size(m,1)
    tmp = crtlFor(m(ii,:)',m(ii+1:end,:)');
    vol(:) = tmp;
    h(ii).fname = ['Othogonalized/' h(ii).fname];
    spm_write_vol(h(ii),vol);
end
