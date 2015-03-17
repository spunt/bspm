% bspm_dartel_norm_func2(input); 
%
% epipat = cell array of patterns for grabbing epi images to norm
% flowfields = flowfields (i.e. u_rc1*) (same length/order as images)
% template = template image (i.e. Template_6.nii)
% voxsize = voxel size for re-sampling (isotropic)
% fwhm = kernel for smoothing (isotropic)
rawdirs             = files('MB*/raw'); 
epipat              = strcat(rawdirs, filesep, 'EP*', filesep, 'bua*nii'); 
flowfields          = files('dartelstuff/MB*/GR*T1*/u_rc1*nii');
input.template      = files('dartelstuff/_templates_/N59*/Template_6*nii');
input.voxsize       = 3; 
input.fwhm          = 6;
allinput = cell(size(flowfields)); 
for i = 1:length(epipat)
    tmp = input; 
    tmp.flowfields = flowfields{i}; 
    tmp.epipat = epipat{i}; 
    allinput{i} = tmp;  
end