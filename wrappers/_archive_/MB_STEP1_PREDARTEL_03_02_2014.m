home; clear all

% paths for relevant folders
path.study = '/Users/bobspunt/Documents/fmri/conte';
qadir = [path.study filesep '_notes_/qa']; mkdir(qadir)

% patterns for finding relevant files/folders 
pattern.subdir = '*CC*';
pattern.epidir = 'EP*';
pattern.t1dir = 'GR*T1*';
pattern.epiimg = 'fad*nii';
pattern.anatimg = 'sad*nii';
fieldpat = 'GR*Field*'; 

% find relevant files folders
subs = files([path.study filesep pattern.subdir]);

% fieldmap options
rot1 = .54*80;
rot2 = .47*32;
rot = [repmat(rot1,length(cellstrfind(subs,'MB_CC')), 1); repmat(rot2,length(cellstrfind(subs,'SB_CC')), 1)];
blip = -1;
method = 'Mark3D';
jacob = 0;
et = [2.55 5.01];

% omit initial volumes
omitpat = {'fad*1-000001*nii' 'fad*2-000002*nii' 'fad*3-000003*nii' 'fad*4-000004*nii'};
bob_spm_omit_vols('MB_CC*/raw/EP*',omitpat);
bob_spm_omit_vols('SB_CC*/raw/EP*',omitpat(1:2));


return
bob_display_message(sprintf('Found %d subjects in %s', length(subs), path.study));
matlabpool 4

parfor s = 1:length(subs)

    subdir = subs{s};
    [ps sub ext] = fileparts(subdir);
    notedir = [subdir filesep 'notes'];

    % grab epi directories
    epidirs = files([subdir filesep 'raw' filesep pattern.epidir]);
    epi_all = cell(epi_all);
    qa_epis = epi_all;
    qa_runnames = epi_all;
    phase_map = epi_all;
    
    for e = 1:length(epidirs)
        
        %% DELETE FIRST N VOLUMES
        omitdir = [epidirs{e} filesep '_omit_'];
        mkdir(omitdir);
        for i = 1:4
            tmp = files([epidirs{e} filesep omitpat{i}]);
            if ~isempty(tmp)
                movefile(char(tmp),omitdir);
            end
        end
        
        %% SLICE TIMING
        epi = files([epidirs{e} filesep pattern.epiimg]);
        hdr = files([epidirs{e} filesep 'dicom*mat']);
        dicominfo = bob_spm_get_dicom_info(char(hdr),0);
        slice_order = dicominfo.sliceinfo.order;
        nslices = dicominfo.sliceinfo.number;
        reference_slice = slice_order(round(nslices/2));
        microtime_resolution = nslices;
        microtime_onset = round(nslices/2);
        TR = dicominfo.parameterinfo.TR/1000;
        bob_spm_slicetime(epi, nslices, TR, slice_order, reference_slice)

        %% FIELDMAP
        fieldmapimg = files([epidirs{e} filesep '*Fieldmap*' filesep pattern.anatimg]);
        mag = fieldmapimg(1:2);
        phase = fieldmapimg(3);
        epi = files([epidirs{e} filesep pattern.epiimg]);
        bob_spm_fieldmap(mag, phase, epi(1), fm, et, rot(s), blip, jacob, method);
        
        %% SAVE SOME VOLUMES
        epi_all{e} = epi;
        phase_map(e) = files([epidirs{e} filesep '*Fieldmap*' filesep 'vdm*nii']);
        [ps name ext] = fileparts(epi{1});
        qa_epis{e} = [ps filesep 'u' name ext];
        [ps name ext] = fileparts(epidirs{e});
        qa_runnames{e} = name;
    
    end    
    
    %% realign and unwarp
    bob_spm_realign_and_unwarp(epi_all, phase_map);
    
    %% compute mean T1
    t1 = files([subdir filesep 'raw' filesep pattern.t1dir filesep pattern.anatimg]);
    meant1dir = [subdir filesep 'raw' filesep 'GR_T1_MEAN'];
    mkdir(meant1dir)
    bob_spm_coregister(t1(1),t1(2));
    meant1 = [meant1dir filesep 'meanT1.nii'];
    bob_spm_imcalc(t1,meant1,'mean');

    %% co-register mean T1 to 
    mean_epi = files([epidirs{1} filesep 'mean*nii']);
    bob_spm_coregister(mean_epi, meant1)

    %% run new segment
    bob_spm_new_segment(meant1);

    %% skullstrip unwarped epis
    for i = 1:length(epidirs)
        epi = files([epidirs{i} filesep 'ua*.nii']);
        bob_spm_bet_epi_batch(epi,.3,'b');
    end

end
matlabpool close

    





