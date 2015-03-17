home; clear all

PRINT = 1;

% paths for relevant folders
path.study = '/Users/bobspunt/Documents/fmri/conte';
path.fmridata = '/Users/bobspunt/Desktop/Dropbox/Bob/Research/Caltech/SURF/tasks/fmri/data';
qadir = [path.study filesep '_notes_/qa']; mkdir(qadir)

% patterns for finding relevant files/folders 
pattern.subdir = '*CC*';
pattern.epidir = 'EP*';
pattern.t1dir = 'GR*T1*';
pattern.epiimg = 'fad*nii';
pattern.anatimg = 'sad*nii';
fieldpat = 'GR*Field*'; 

% epi initial volumes to omit
omitpat = {'fad*1-000001*nii' 'fad*2-000002*nii' 'fad*3-000003*nii' 'fad*4-000004*nii'};

% fieldmap options
fm(1).blip = -1;
fm(1).method = 'Mark3D';
fm(1).jacob = 0;
fm(1).et = [2.55 5.01];
fm(1).rot = .47*32;
fm(2).blip = -1;
fm(2).method = 'Mark3D';
fm(2).jacob = 0;
fm(2).et = [2.55 5.01];
fm(2).rot = .54*80;

% slice time stuff
% header information
dicominfo = bob_spm_get_dicom_info('LOI/dicom_header.mat');
slice_order = dicominfo.sliceinfo.order;
nslices = dicominfo.slic         der(round(nslices/2));
microtime_resolution = nslices;
microtime_onset = round(nsl ices/2);  
TR = dicominfo.parameterinfo.TR/1000;
epi = files('LOI/f*nii');
bob_spm_slicetime(epi, nslices, TR, slice_order, reference_slice)



% find relevant files folders
subs = files([path.study filesep pattern.subdir]);
bob_display_message(sprintf('Found %d subjects in %s', length(subs), path.study));

matlabpool 4
parfor s = 1:length(subs)

    subdir = subs{s};
    [ps sub ext] = fileparts(subdir);
    notedir = [subdir filesep 'notes'];
    bob_display_message(sprintf('Working on: %s', subs{s}), '-');

    
    

    % grab epi directories
    epidirs = files([subdir filesep 'raw' filesep pattern.epidir]);
    fieldmapdirs = files([subdir filesep 'raw' filesep fieldpat]);
    mag = files([fieldmapdirs{1} filesep '*.nii']);
    phase = files([fieldmapdirs{2} filesep '*.nii']);
    mag_phase = files([subdir filesep 'raw' filesep fieldpat filesep '*.nii']);

    nrun = length(epidirs);
    epi_all = cell(nrun,1);
    epi_first = cell(nrun,1);
    qa_epis = cell(nrun,1);
    qa_runnames = cell(nrun,1);
    for i = 1:length(epidirs)

        % delete initial volumes
        omitdir = [epidirs{i} filesep '_omit_'];
        mkdir(omitdir);
        for o = 1:length(omitpat)
            tmp = files([epidirs{i} filesep omitpat{o}]);
            if ~isempty(tmp)
                movefile(char(tmp),omitdir);
            end
        end
        tmp = files([epidirs{i} filesep '*.nii']);
        epi_all{i} = tmp;
        epi_first{i} = tmp{1};
        [path name ext] = fileparts(tmp{1});
        qa_epis{i} = [path filesep 'u' name ext];
        [path name ext] = fileparts(epidirs{i});
        qa_runnames{i} = name;
    end

    % compute phase map (vdm) - bob_spm_fieldmap(magnitude_image, phase_image, epi_images, echo_times, total_epi_readout_time, blip, jacobTAG, method)
    bob_spm_fieldmap(mag, phase, epi_first, et, rot, blip, jacob, method);

    % run realign and unwarp - bob_spm_realign_and_unwarp(epi_images, phase_map)
    if length(epidirs)>1
        phase_map = files([fieldmapdirs{2} filesep 'vdm*run*nii']);
    else
        phase_map = files([fieldmapdirs{2} filesep 'vdm*nii']);
    end
    bob_spm_realign_and_unwarp(epi_all, phase_map)

    % compute mean T1
    t1 = files([subdir filesep 'raw' filesep pattern.t1dir filesep pattern.anatimg]);
    meant1dir = [subdir filesep 'raw' filesep 'GR_T1_MEAN'];
    mkdir(meant1dir)
    bob_spm_coregister(t1(1),t1(2));
    meant1 = [meant1dir filesep 'meanT1.nii'];
    bob_spm_imcalc(t1,meant1,'mean');

    % co-register mean T1 to 
    mean_epi = files([epidirs{1} filesep 'mean*nii']);
    bob_spm_coregister(mean_epi, meant1)

    % run new segment
    bob_spm_new_segment(meant1);

    % skullstrip unwarped epis
    for i = 1:length(epidirs)
        epi = files([epidirs{i} filesep 'uf*.nii']);
        bob_spm_bet_epi_batch(epi,.3,'b');
    end

    % print images for evaluation
    if PRINT
        images = [t1; qa_epis'];
        captions = ['T1'; qa_runnames'];
        xpos = 0;
        for p = 1:length(xpos)
            bob_spm_checkreg(images, captions, [xpos(p) 40 0]);
            saveas(gcf, sprintf('%s/%s_REG_%d.jpg', qadir, sub, p), 'jpg');
        end
    end

end
matlabpool close

    





