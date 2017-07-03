function bramila_glm1st(subject,dataroot,outdirname)
addpath(genpath('/m/nbe/scratch/braindata/shared/TouchHyperScan/spm12/spm12b'));
addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila'));
% load cfg and unpack parameters
load(sprintf('%s/%s/%s.mat',dataroot,subject,outdirname),'cfg');
TR = cfg.TR;
sesses = cfg.sesses;
hpf = cfg.hpf;
incmoves = cfg.incmoves;
inputfile = cfg.inputfile;
motionfile = cfg.motionfile;
cnames = cfg.cnames;
regressor_unit = cfg.regressor_unit;
onsets = cfg.onsets;
durations = cfg.durations;

if ~isfield(cfg,'AR'), cfg.AR = 'none'; end
AR = cfg.AR;

% where would analysis result go
outputdir = fullfile(dataroot,subject,outdirname);
mkdir(outputdir);
%number of blocks
nsess = length(sesses);
%condition names
ncond = length(cnames);
% folder where the fun happens
for sess = 1:nsess
    sessdata = fullfile(dataroot,subject,sesses{sess});
%     %% Split 4D volume into 3D volumes. When file is larger than 2.1gb, you have to split it
%     % assuming sessdata is folder with preprocessed data for given run/session
%     % assuming inputfile is epi.nii
%     if exist(fullfile(sessdata,'split'),'dir') ~= 7 || length(dir(fullfile(sessdata,'split')))==2% if folder does not exist or is empty
%         mkdir(fullfile(sessdata,'split'));
%         spm_file_split(fullfile(sessdata,inputfile),fullfile(sessdata,'split'));
%     end
%     %% List files
%     clear files
%     clear ffiles
%     files = spm_select('List',fullfile(sessdata,'split'),'^*\.nii$');
%     for f = 1:length(files)
%         ffiles{f,1} = [fullfile(sessdata,'split',files(f,:)),',1'];
%     end
%     % Throw them into the batch
%     matlabbatch{1}.spm.stats.fmri_spec.sess(sess).scans = ffiles;
    %% SPM12 seems to work just fine
    [ffiles,~] = spm_select('ExtFPList',sessdata,inputfile,Inf);
    % Throw them into the batch
    matlabbatch{1}.spm.stats.fmri_spec.sess(sess).scans = cellstr(ffiles);        
    %% Regressors
    % Create and throw into the batch
    for c = 1:ncond
        matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).name = cnames{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).onset = onsets{sess}{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).duration = durations{sess}{c};
        % no temporal or parametric modulations here
        matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).tmod = 0;
        if cfg.pmod == 1
            for p = 1:length(cfg.pmodname)
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).pmod(p).name = cfg.pmodname{p};
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).pmod(p).param = cfg.pmodparam{sess}{c}{p};
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).pmod(p).poly = 1;
            end
        end
    end
    %% Nuisance
    matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(sess).regress = struct('name', {}, 'val', {});
    if incmoves==1
        R = load([sessdata '/' motionfile]);
        if length(R) ~= length(ffiles)
            disp('WARNING: the length of motion parameters is not equal to length of files, need to trim')
        end
        save([sessdata '/' motionfile '.mat'],'R');
        matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi_reg = {[sessdata '/' motionfile '.mat']};        
    else
        matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi_reg = {''};
    end
    
    %% High-pass filter
    matlabbatch{1}.spm.stats.fmri_spec.sess(sess).hpf = hpf;
end
%% Put together the other parts of batch
% where would results be
matlabbatch{1}.spm.stats.fmri_spec.dir = {outputdir};
% seconds or scans
matlabbatch{1}.spm.stats.fmri_spec.timing.units = regressor_unit;
% TR
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
% microtiming (ignore)
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
% factorial model (ignore)
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
% derivatives of hrf 
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;
matlabbatch{1}.spm.stats.fmri_spec.mask = {'/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/external/MNI152_T1_2mm_brain_mask.nii,1'};
matlabbatch{1}.spm.stats.fmri_spec.cvi = AR;
%% Model estimation barch
matlabbatch{2}.spm.stats.fmri_est.spmmat = {[outputdir '/SPM.mat']};
% Classical or Bayesian? Dunno...
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
%% Contrast manager (can be ignored in future, if you want to understand it better)
if ~isempty(cfg.predefined_contrast)
    % If we defined contrasts in advance
    for cond = 1:length(cfg.predefined_contrast)
        contrst{cond} = cfg.predefined_contrast{cond};
        contrstname{cond} = cfg.predefined_contrast_name{cond};        
    end
else
    if ncond > 1
        contrst = cell(ncond*2,1);
        contrstname = contrst;
        % main effect is easy to manage
        for cond = 1:ncond
            % build the base of the contrast    
            contrst{cond} = zeros(1,ncond);
            contrst{cond}(cond) = 1;
            contrstname{cond} = [cnames{cond} ' Main Effect'];
        end
        % contrast (one vs all) is more sophisticated
        % cond + ncond to not overwrite anything
        for cond = 1:ncond
            % build the base of the contrast
            consbase = ones(1,ncond);
            consbase = consbase*-1;
            consbase(cond) = ncond - 1;        
            contrst{cond+ncond} = consbase;
            contrstname{cond+ncond} = [cnames{cond} ' vs All'];
        end
    elseif ncond == 1
        contrst = cell(1,1);
        contrstname = contrst;
        % main effect is easy to manage
        contrst{1} = 1;
        contrstname{1} = [cnames{1} ' Main Effect'];
    else
        disp('What is wrong with your cnames?')
    end
end
matlabbatch{3}.spm.stats.con.spmmat = {[outputdir '/SPM.mat']};
matlabbatch{3}.spm.stats.con.delete = 0;
for ct = 1:length(contrst)
    matlabbatch{3}.spm.stats.con.consess{ct}.tcon.name = contrstname{ct};
    matlabbatch{3}.spm.stats.con.consess{ct}.tcon.convec = contrst{ct};
    matlabbatch{3}.spm.stats.con.consess{ct}.tcon.sessrep = 'repl';
end
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
