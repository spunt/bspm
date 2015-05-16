% ======================================
% a wrapper for setting up an fMRI subject using dicom data
% ======================================

% do flags
do_convert = 0;
do_only_convert = 0;
do_print = 0;
do_realign = 0;
do_slice_timing = 0;
do_norm = 0;

% paths and patterns for relevant files and folders
studydir = '/Users/bobspunt/Documents/fmri/lois';
qadir = [studydir filesep 'qa']; mkdir(qadir)
subpat = 'RA0548*';
epipat = 'EP*';
t1pat = 'GR*T1*';
fieldpat = 'GR*Field*'; 

% fieldmap options
blip = -1;
method = 'Mark3D';
jacob = 0;          % flag for running jacobian modulation
et = [2.55 5.01];  % echo times (ms)
rot = .47*32;       % epi readout time (ms)

% bad scan options
skip = 2; % number of initial volumes defined as bad
thresh.globalsignalchange = 3; % in standard deviations
thresh.motionmag = .5; % in mm (trans) or degrees (rot)
rpTAG = 1; % 1 = include, 0 = omit;
filename = 'badscan_skip2_gs300.txt';

% normalization/smoothing options
voxsize = 3;
fwhm = 8;

% DO IT
% ----------------

subs = files([studydir filesep subpat]);
bspm_display_message(sprintf('Setting up %d subjects in %s', length(subs), studydir));

% for s = 1:length(subs)
for s = 8
    
    subdir = subs{s};
    [path sub ext] = fileparts(subdir);
    notedir = [subdir filesep 'notes'];
    bspm_display_message(sprintf('Working on: %s', subs{s}), '-');

    if do_convert
        % convert dicom to nii - bspm_convert_dcm(subDIR, format)
        bspm_convert_dcm(subdir, 'nii');
    end
    if do_only_convert
        continue
    end
    
    % grab important images
    fieldmapdirs = files([subdir filesep 'raw' filesep fieldpat]);
    mag = files([fieldmapdirs{1} filesep '*.nii']);
    phase = files([fieldmapdirs{2} filesep '*.nii']);
    mag_phase = files([subdir filesep 'raw' filesep fieldpat filesep '*.nii']);
    epidirs = files([subdir filesep 'raw' filesep epipat]);
    epi_all = {};
    epi_first = {};
    qa_epis = {};
    qa_runnames = {};
    for i = 1:length(epidirs)
        tmp = files([epidirs{i} filesep '*.nii']);
        epi_all{i} = tmp;
        epi_first{i} = tmp{1};
        [path name ext] = fileparts(tmp{1});
        qa_epis{i} = [path filesep 'u' name ext];
        [path name ext] = fileparts(epidirs{i});
        qa_runnames{i} = name;
    end
    t1 = files([subdir filesep 'raw' filesep t1pat filesep '*.nii']);    

    if do_realign
        % compute phase map (vdm) - bspm_fieldmap(magnitude_image, phase_image, epi_images, echo_times, total_epi_readout_time, blip, jacobTAG, method)
        bspm_fieldmap(mag, phase, epi_first, et, rot, blip, jacob, method);

        % run realign and unwarp - bspm_realign_and_unwarp(epi_images, phase_map)
        phase_map = files([fieldmapdirs{2} filesep 'vdm*nii']);
        bspm_realign_and_unwarp(epi_all, phase_map)
    end
    
    % co-register anatomial to skull-stripped mean epi
    mean_epi = files([epidirs{1} filesep 'mean*nii']);
%     mean_epi_brain = [epidirs{1} filesep 'mean_epi_brain'];
%     command = sprintf('bet %s %s -f 0.3', char(mean_epi), mean_epi_brain); system(command);
%     gunzip([mean_epi_brain '.nii.gz']);
%     bspm_coregister(files([epidirs{1} filesep 'mean*brain*nii']), t1)
    bspm_coregister(mean_epi, t1)
    
    % print images for evaluation
    if do_print
        images = [t1; qa_epis'];
        captions = ['T1'; qa_runnames'];
        xpos = 0;
        for p = 1:length(xpos)
            bspm_checkreg(images, captions, [xpos(p) 40 0]);
            saveas(gcf, sprintf('%s/%s_REG_%d.jpg', qadir, sub, p), 'jpg');
        end
    end
        
    % run badscan, then slice time correct if desired
    for i = 1:length(epidirs)
        
        bspm_badscan(epidirs{i}, 'u*nii', skip, thresh, rpTAG, filename)
        
        if do_slice_timing
            % get info for slice timing
            load([epidirs{i} filesep 'dicom_header.mat']);
            slicetimes = hdr{1}.Private_0019_1029';
            TR = hdr{1}.RepetitionTime/1000;
            nslices = length(slicetimes);
            slicetimes(:,2) = 1:nslices;
            slicetimes = sortrows(slicetimes,1);
            slice_order = slicetimes(:,2);
            reference_slice = slicetimes(round(nslices/2),2);
            epi_images = files([epidirs{i} filesep 'u*nii']);
            bspm_slicetime(epi_images, nslices, TR, slice_order, reference_slice)
            epi_prefix = 'au';
        else
            epi_prefix = 'u';
        end

    end

    % grab all epis (across sessions) for normalization
    allepi{s} = files([subdir filesep 'raw' filesep epipat filesep epi_prefix '*.nii']);

end

if do_norm

    % run new segment
    allt1 = files([studydir filesep subpat filesep 'raw' filesep t1pat filesep '*.nii']);    
    bspm_new_segment(allt1);

    % create dartel template
    anatpat = [studydir filesep subpat filesep 'raw' filesep t1pat];
    bspm_dartel_create_template(anatpat, 'nii');

    % normalize anatomicals and functionals
    bspm_dartel_normalize(allt1, anatpat, 'nii', 1, 0);
    bspm_dartel_normalize(allepi, anatpat, 'nii', voxsize, fwhm);

end





 
 
 
 
 
 
 
 
 
 
 
 
