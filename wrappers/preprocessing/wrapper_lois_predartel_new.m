home; clear all

% | options
runit.omitvols      = 0;
runit.slicetime     = 1;
opt.slice_times     = 0; % 1 do slice-timing using actual times, 0 will do using order
opt.coreg_epi2t1    = 0; % 0 will coreg t1 to mean EPI; 1 will coreg all EPI to t1
runit.segment       = 0;

% | paths for relevant folders
path.study  = '/Users/bobspunt/Dropbox/Writing/Empirical/ASD/data';
if regexp(pwd,'/home/spunt'), path.study = regexprep(path.study, '/Users/bobspunt/Documents/fmri', '/home/spunt/data'); end
path.qa     = fullfile(path.study, '_qa_');
if ~exist(path.qa, 'dir'), mkdir(path.qa); end

% | patterns for finding relevant files/folders (relative to subject dir)
pattern = struct( ...
    'subdir',   'RA0979_LOIS_A_PEERS*',          ...
    'epidir',   'EP*',          ...
    't1dir',    'GR*T1*',       ...
    'fmdir',    'GR*Field*',    ...
    'anatimg',  'sa*nii',       ...
    'epiimg',   'fa*nii',       ...
    'refdcm',   '*lois*dcm'     ...
    );

% | Field Map Parameters
rot     = .47*32;   % SURF: .54*80
blip    = -1;
method  = 'Mark3D';
jacob   = 0;
ets     = [2.55 5.01];  % SURF: [2.15 4.61]

% | reference dicom image
refdcm      = bspm_dicomref(pattern.refdcm);

% | relevant directories
subdirs      = files(fullfile(path.study, pattern.subdir));
omitpat      = [];
if ~isempty(omitpat), subdirs(cellstrfind(subdirs, omitpat)) = []; end

% | Omit Initial Volumes
if runit.omitvols
    omitpat = {'fad*_00001.nii' 'fad*_00002.nii'};
    bspm_omit_vols(fullfile(path.study, pattern.subdir, 'raw', pattern.epidir), omitpat);
end

% | Build Jobs Subject-by-Subject
allsub = cell(length(subdirs), 1);
count   = 0;
for s = 1:length(subdirs)

    subdir          = subdirs{s};
    [ps, sub, ext]  = fileparts(subdir);

    % | grab epi directories
    epidirs     = files([subdir filesep 'raw' filesep pattern.epidir]);
    epi_all     = cell(size(epidirs));
    qa_epis     = epi_all;
    qa_runnames = epi_all;
    phase_map   = epi_all;
    uaepi       = [];
    for e = 1:length(epidirs)

        % | Define EPIs
        epi = files([epidirs{e} filesep pattern.epiimg]);
        epi = bspm_check_filenames(epi);
        [epip, epin, epie] = cellfun(@fileparts, epi, 'unif', false);
        pat1 = {'' 'a'};
        pat2 = {'u' 'ua'};
        epi_st = strcat(epip, filesep, pat1{runit.slicetime+1}, epin, epie);
        epi_uw = strcat(epip, filesep, pat2{runit.slicetime+1}, epin, epie);

        % | Slice Timing
        if runit.slicetime
            count       = count + 1;
            matlabbatch(count)   = bspm_slicetime(epi, refdcm, opt.slice_times+1);
        end

        % | Field Map
        fieldmapimg = files([epidirs{e} filesep '*Fieldmap*' filesep pattern.anatimg]);
        mag         = fieldmapimg(1:2);
        phase       = fieldmapimg(3);
        count       = count + 1;
        matlabbatch(count)   = bspm_fieldmap(mag, phase, epi_st(1), rot, 'ets', ets, 'blip', blip, 'jacob', jacob, 'method', method);

        % | Save Volumes for Unwarp
        epi_all{e} = epi_st;
        [pmapp, pmapn, pmape] = fileparts(phase{1});
        phase_map{e} = strcat(pmapp, filesep, 'vdm5_sc', pmapn, pmape);
        qa_epis{e} = epi_uw{1};
        if e==1
            [ep, en, ee] = fileparts(epi_uw{1});
            mean_epi = strcat(ep, filesep, 'mean', en, ee);
        end
        [ps name ext] = fileparts(epidirs{e});
        qa_runnames{e} = name;

        % | Save FileNames for Coregistration to T1
        uaepi = [uaepi; epi_uw];

    end

    % | Realign and Unwarp
    count = count + 1;
    matlabbatch(count) = bspm_realign_and_unwarp(epi_all, phase_map);

    % | Co-register
    t1 = files([subdir filesep 'raw' filesep pattern.t1dir filesep pattern.anatimg]);
    count = count + 1;
    if opt.coreg_epi2t1
        % | EPIs to T1
        matlabbatch(count) = bspm_coregister(t1, mean_epi, uaepi);
    else
        % | T1 to Mean EPI
        matlabbatch(count) = bspm_coregister(mean_epi, t1);
    end

    % | Segment T1
    if runit.segment
        count = count + 1;
        matlabbatch(count) = bspm_segment(t1);
    end

end

% | Save Job
jobname = fullfile(pwd, sprintf('job_pp_%s.mat', strtrim(datestr(now,'mmm_DD_YYYY'))));
save(jobname, 'matlabbatch');

return

% | Add to Existing Template
rc1 = files([subdir filesep 'raw' filesep pattern.t1dir filesep 'rc1*nii']);
rc2 = files([subdir filesep 'raw' filesep pattern.t1dir filesep 'rc2*nii']);
template = = files(fullfile(path.study, '**', 'Template*6.ni*'));
templatedir = fileparts(template);
pigz(fullfile(templatedir, '*gz'));
matlabbatch = bspm_dartel_existing_template(rc1, rc2, templatedir)




