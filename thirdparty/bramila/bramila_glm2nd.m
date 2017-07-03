% 2nd level (random effects) modeling in SPM
clear all
close all
addpath(genpath('/m/nbe/scratch/braindata/shared/TouchHyperScan/spm12/spm12b'));
dataroot = '/m/nbe/scratch/braindata/shared/GraspHyperScan';
% outdirname is the folder where 1st level results were saved (see bramila_glminit.m)
outdirname = {
'1stlevel';
};
% Where to save the 2nd level models
resultsfolder = 'Observers';
% Subject list
subjects = {
'Fanny_Observer';
'Fanny_Observer_1';
'Fanny_Observer_2';
'Fanny_Observer_3';
'Sonya_Observer';
'Sonya_Observer_1';
'Sonya_Observer_2';
'Sonya_Observer_3';
'Sonya_Observer_4';
'Sonya_Observer_5';
'Sonya_Observer_6';
'Sonya_Observer_7';
'Sonya_Observer_8';
'Sonya_Observer_9';
'Sonya_Observer_10';
'Sonya_Observer_11'
}; 
%% model specification, take same contrasts as in 1st level
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); % if you have some covariates, add one per subject
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/external/MNI152_T1_2mm_brain_mask.nii,1'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
% Cell array of cells
% In case you have several analyses
% In case you have one, just have a cell within cell, which has the
% contrasts you ran in 1st level model
contrastnames = {
{'Ball','Pen','Slap','Point','Ball vs All','Pen vs All','Slap vs All','Point vs All'};
};
% how files representing each contrast will be called in 1st lvl
contrastnr = {'con_0001.nii','con_0002.nii','con_0003.nii','con_0004.nii','con_0005.nii','con_0006.nii','con_0007.nii','con_0008.nii'};
%% Execute
for a = 1:length(contrastnames)
    contrasts = contrastnames{a};   
    for cons = 1:length(contrasts);
        % specify analysis directory
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {sprintf('%s/%s/%s/SPM.mat',dataroot,resultsfolder,contrasts{cons})};
        anadir = fullfile(dataroot,resultsfolder,contrasts{cons});
        if ~exist(anadir), mkdir(anadir); end
        matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(anadir);
        % specify files from 1st lvl
        clear files
        for s = 1:length(subjects)
            funcdir = fullfile(dataroot,subjects{s},outdirname{a});
            files{s,1} = fullfile(funcdir,contrastnr{cons}); 
        end
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = files;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = contrasts{cons};
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = 1;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
        matlabbatch{3}.spm.stats.con.delete = 0;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.spmmat = {sprintf('%s/%s/%s/SPM.mat',dataroot,resultsfolder,contrasts{cons})};
        %Result report, partly tested
        %    matlabbatch{4}.spm.stats.results.spmmat = {sprintf('%s/%s/%s/SPM.mat',dataroot,resultsfolder,contrasts{cons})};
        %    %Contrast number
        %    matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
        %    %FDR, FWE or none
        %    matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
        %    %pvalue
        %    matlabbatch{4}.spm.stats.results.conspec.thresh = 0.0500;
        %    %cluster size
        %    matlabbatch{4}.spm.stats.results.conspec.extent = 0;
        %    matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
        %    matlabbatch{4}.spm.stats.results.print = 'pdf';
        %    %Thresholded nifti output file basename
        %    matlabbatch{4}.spm.stats.results.write.tspm.basename = 'thresholded_SPM';
        if exist([anadir '/SPM.mat'],'file') ~= 0
            delete([anadir '/SPM.mat'])
        end
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);    
    end
end



