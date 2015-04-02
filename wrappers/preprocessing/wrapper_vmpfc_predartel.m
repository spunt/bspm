home; clear all

% | options
runit.omitvols      = 1; 
runit.slicetime     = 1;
opt.slice_times     = 0; % 1 do slice-timing using actual times, 0 will do using order
opt.coreg_epi2t1    = 0; % 0 will coreg t1 to mean EPI; 1 will coreg all EPI to t1
runit.segment       = 0;

% | paths for relevant folders
path.study  = '/Users/bobspunt/Documents/fmri/vmpfc';
path.qa     = fullfile(path.study, '_notes_', 'qa');
if ~exist(path.qa, 'dir'), mkdir(path.qa); end

% | patterns for finding relevant files/folders (relative to subject dir)
pattern = struct(                   ...
    'subdir',   'sub0747*',             ...
    'epidir',   'EP*LOI_2*',              ...
    't1dir',    'GR*T1*',           ...
    'fmdir',    'GR*Field*LOI*',        ...
    'anatimg',  's*nii',                ...
    'epiimg',   'f*nii',                ...
    'refdcm',   fullfile(fileparts(mfilename('fullpath')), 'ref*dicom**vmpfc.dcm')     ...
    );

% | reference dicom image
refdcm      = files(pattern.refdcm);
dcm         = bspm_get_dicom_info(refdcm); 
nlines      = dcm.parameterinfo.matrixsize(1);
ipatfactor  = dcm.sequenceinfo.ipatfactor; 

% | Field Map Parameters
rot     = .47*(nlines/ipatfactor);  % echo spacing * # of lines of data acquired
blip    = -1;
method  = 'Mark3D';
jacob   = 0;
ets     = [2.55 5.01]; 

% | relevant directories
subdirs      = files(fullfile(path.study, pattern.subdir)); 
omitpat      = [];
if ~isempty(omitpat), subdirs(cellstrfind(subdirs, omitpat)) = []; end

% | Omit Initial Volumes
if runit.omitvols
    omitpat = {'f*00001-000001*nii' 'f*00002-000002*nii'};
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

        % | Save Volumes for Unwarp
        epi_all{e} = epi_st; 
%         [pmapp, pmapn, pmape] = fileparts(phase_map{1}); 
%         phase_map{e} = strcat(pmapp, filesep, 'vdm5_sc', pmapn, pmape); 
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
    
    
    
    % | Field Map 
%     epi_st      = [epi_all{1}(1); epi_all{2}(1)]
    epi_st      = epi_all{1}(1); 
    fieldmapimg = files([subdir filesep 'raw' filesep pattern.fmdir filesep pattern.anatimg]);
    mag         = fieldmapimg(1:2);
    phase       = fieldmapimg(3);
    count       = count + 1;
    matlabbatch(count)   = bspm_fieldmap(mag, phase, epi_st, rot, 'ets', ets, 'blip', blip, 'jacob', jacob, 'method', method);

    return
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





