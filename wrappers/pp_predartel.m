home; clear all

% | paths for relevant folders
path.study  = '/Users/bobspunt/Documents/fmri/lois';
path.qa     = fullfile(path.study, '_notes_', 'qa');
if ~exist(path.qa, 'dir'), mkdir(path.qa); end

% | patterns for finding relevant files/folders (relative to subject dir)
pattern = struct( ...
    'subdir',   'RA*',          ...
    'epidir',   'EP*',          ...
    't1dir',    'GR*T1*',       ...
    'fmdir',    'GR*Field*',    ...
    'anatimg',  'sa*nii',       ...
    'epiimg',   'fa*nii'        ...
    );

% | relevant directories
subdirs      = files(fullfile(path.study, pattern.subdir)); 
omitpat      = 'RA0082'; 
subdirs(cellstrfind(subdirs, omitpat)) = []; 

% | reference dicom image
refdcm      = files('dicom_ref*.dcm');

% | fieldmap options
rot     = .47*32;
blip    = -1;
method  = 'Mark3D';
jacob   = 0;
ets     = [2.55 5.01];

% | omit initial volumes
omitpat = {'fad*1-000001*nii' 'fad*2-000002*nii'};
bspm_omit_vols(fullfile(path.study, pattern.subdir, 'raw', pattern.epidir), omitpat); 

% | Build Jobs
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

    for e = 1:length(epidirs)
                
        % | Define EPIs
        epi = files([epidirs{e} filesep pattern.epiimg]);
        [epip, epin, epie] = cellfun(@fileparts, epi, 'unif', false); 
        epi_st = strcat(epip, filesep, 'a', epin, epie);
        epi_uw = strcat(epip, filesep, 'ua', epin, epie); 

        % | Slice Timing
        count       = count + 1;
        mb(count)   = bspm_slicetime(epi, refdcm, 1); 

        % | Field Map 
        fieldmapimg = files([epidirs{e} filesep '*Fieldmap*' filesep pattern.anatimg]);
        mag         = fieldmapimg(1:2);
        phase       = fieldmapimg(3);
        count       = count + 1;
        mb(count)   = bspm_fieldmap(mag, phase, epi_st(1), rot, 'ets', ets, 'blip', blip, 'jacob', jacob, 'method', method);
        
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
    
    end    
    
    % | Realign and Unwarp
    count = count + 1; 
    mb(count) = bspm_realign_and_unwarp(epi_all, phase_map); 
    
    % | Co-register T1 to Mean EPI
    t1 = files([subdir filesep 'raw' filesep pattern.t1dir filesep pattern.anatimg]);
    count = count + 1; 
    mb(count) = bspm_coregister(mean_epi, t1); 

    % | Segment T1
    count = count + 1; 
    mb(count) = bspm_segment(t1);
    
end

% | Save Job
jobname = fullfile(pwd, sprintf('job_pp_%s.mat', strtrim(datestr(now,'mmm_DD_YYYY')))); 
save(jobname, 'mb'); 





