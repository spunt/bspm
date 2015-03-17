clear all; home; 

path.study      = '/home/spunt/data/lois';
pattern.subdir = 'MB_RA0*';
pattern.epidir = 'EP*LOI*';
pattern.t1dir = 'GR*T1*';
pattern.fielddir = 'GR*Field*';
pattern.epiimg = 'fad*nii';
pattern.epiimg2 = 'afad*nii';
pattern.anatimg = 'sad*nii';

%% define pattern inputs
subdirs = files([path.study filesep pattern.subdir]);
rawpat = [path.study filesep pattern.subdir filesep 'raw'];
epidirs = files([rawpat filesep pattern.epidir]);
fielddirs = files([rawpat filesep pattern.epidir filesep pattern.fielddir]);
t1dirs = files([rawpat filesep pattern.t1dir]);
t1pat = [rawpat filesep pattern.t1dir];
epipat_raw = strcat(epidirs,[filesep pattern.epiimg]);
epipat_slicetimed = strcat(epidirs,[filesep 'a' pattern.epiimg]);
epipat_unwarped = strcat(epidirs,[filesep 'ua' pattern.epiimg]);
omitpat = {'fad*1-000001*nii' 'fad*2-000002*nii' 'fad*3-000003*nii' 'fad*4-000004*nii'};
folderpat = [rawpat filesep pattern.epidir];
input = bob_cell2struct([{omitpat} {folderpat}],{'omitpat' 'folderpat'});
bob_spm_omit_vols2(input)

%% SLICE TIMING
hdr = files([epidirs{1} filesep 'dicom*mat']);
dicominfo = bob_spm_get_dicom_info(char(hdr),0);
slice_order = dicominfo.sliceinfo.order;
nslices = dicominfo.sliceinfo.number;
reference_slice = slice_order(round(nslices/2));
microtime_resolution = nslices;
microtime_onset = round(nslices/2);
TR = dicominfo.parameterinfo.TR/1000;
input.slice_order = slice_order;
input.nslices = nslices;
input.reference_slice = reference_slice;
input.TR = TR;
allinput.slicetime = cell(length(epidirs),1);
for e = 1:length(epidirs)
    input.epipat = epipat_raw{e};
    allinput.slicetime{e} = input;
end
clear input

%% FIELDMAP
input.blip = -1;
input.method = 'Mark3D';
input.jacob = 0;                % flag for running jacobian modulation
input.echotimes = [2.55 5.01];  % echo times (ms)
input.epirot = .54*80;          % epi readout time (ms)
allinput.fieldmap = cell(length(epidirs),1);
for e = 1:length(epidirs)
    input.epipat = epipat_slicetimed{e};
    fmimg = files([epidirs{e} filesep pattern.fielddir filesep pattern.anatimg]);
    input.magimg = fmimg(1:2);
    input.phaseimg = fmimg(3);
    allinput.fieldmap{e} = input;
end
clear input

%% REALIGN & UNWARP
allinput.unwarp = cell(length(subdirs),1);
for s = 1:length(subdirs)
    input.phasepat = [];
    input.epipat = [];
    subepidirs = files([subdirs{s} filesep 'raw' filesep pattern.epidir]);
    subepipat = strcat(subepidirs,[filesep 'a' pattern.epiimg]);
    subfieldmapdirs = files([subdirs{s} filesep 'raw' filesep pattern.epidir filesep pattern.fielddir]);
    input.phasepat = strcat(subfieldmapdirs,[filesep 'vdm*nii']);
    input.epipat = subepipat;
    allinput.unwarp{s} = input;
end
clear input

allinput.coreg = cell(length(subdirs),1);
allinput.segment = cell(length(subdirs),1);
for s = 1:length(subdirs)

    input.reference = files([subdirs{s} filesep 'raw' filesep pattern.epidir filesep 'mean*nii']);
    input.source = files([subdirs{s} filesep 'raw' filesep pattern.t1dir filesep pattern.anatimg]);
    input2.img = input.source;
    allinput.coreg{s} = input;
    allinput.segment{s} = input2;
end
clear input
clear input2

return

% bob_submit2cluster(@bob_spm_slicetime2, allinput.slicetime, 30);
% bob_submit2cluster(@bob_spm_fieldmap2, allinput.fieldmap, 30);
% bob_submit2cluster(@bob_spm_realign_and_unwarp2, allinput.unwarp, 30);

% bob_submit2cluster(@bob_spm_coregister2,allinput.coreg,30);
% allinput.segment(20) = [];
% bob_submit2cluster(@bob_spm_new_segment2,allinput.segment,30);
% anatpat = {t1pat};
% bob_submit2cluster(@bob_spm_dartel_create_template,{anatpat},30);


































return


flowfields = files([t1pat filesep 'u_rc*nii']);
voxsize = 2;
fwhm = 0;
allinput = cell(length(subdirs),1);
template = files([path.study filesep '_tem*/T*6.nii']);
for e = 1:length(subdirs)
    input.epipat = [subdirs{e} filesep 'raw' filesep pattern.epidir filesep 'ua*nii'];
%     input.thresh.globalsignalchange = 2.5;
%     input.thresh.motionmag = .5;
    
    input.flowfields = flowfields(e);
    input.template = template;
    input.voxsize = voxsize; 
    input.fwhm = fwhm;
    allinput{e} = input;
end
input = allinput; 


