function bramila_preprocessor(cfg)
% BRAMILA_PREPROCESSOR - executes the single subject preprocessing pipeline.
%   - Usage:
%   bramila_preprocessor(cfg,gridtype)
%   - Input:
%   cfg is a struct with following parameters
%		### ENRICO TO ADD MORE ###
%       cfg.subj = cfg is a struct containing all the parameters from 
%       BRAMILA_PREPRO_DEMO.m
%   gridtype = boolean, 1 if slurm, 0 if else.
%   - Notes:
%   Refer to wiki page:
%   https://wiki.aalto.fi/pages/viewpage.action?pageId=92445275 for
%   detailed documentation
%	Last edit: EG 2014-07-14


% adding some toolboxes that we might need
addpath(('/m/nbe/scratch/braindata/shared/toolboxes/spm12b/'));
addpath(('/m/nbe/scratch/braindata/shared/toolboxes/afni_matlab/'));



% Set fsldir
if ischar(cfg) % in case we use triton, cfg will be a path to cfg
    load(cfg);
    setenv('FSLDIR','/share/apps/fsl/5.0.9/fsl/') % it was setenv('FSLDIR','/share/apps/fsl/4.1.9');

else
    % First, some library problems
    [~,hname] = system('hostname');
	if isempty(hname)
		% sometimes (???) hname is empty... no idea why
		hname='whatever';
	end
	        disp(['You are on machine ' hname])

    hname = hname(1:4); % to remove the nasty space in the end
    
	if strcmp(hname,'fn01') == 1 || strcmp(hname,'joun') == 1 % if fn01 - use triton FSLDIRa %%% THIS NEEDS TO BE UPDATED FOR INTERACTIVE SESSIONS
        setenv('FSLDIR','/share/apps/fsl/5.0.9/fsl/') % it was setenv('FSLDIR','/share/apps/fsl/4.1.9');
    else
	    setenv('FSLDIR','/usr/share/fsl/5.0');
		% the below c++ libraries don't exist anymore on dione... 
    %    setenv('LD_PRELOAD','/usr/lib/gcc/x86_64-linux-gnu/4.6.3/libstdc++.so'); % some library that my matlab couldnt find on local machine, but can find on triton
    end
	% this was missing, I am not sure how it was possible to work before
	system('source $FSLDIR/etc/fslconf/fsl.sh')
end

if(~isfield(cfg,'overwrite'))
	cfg.overwrite = 0;
end

% the below was hardcoded string. It could check for existence of the path
addpath(genpath(cfg.bramilapath));
% outputs the git release used
githashfile=[cfg.bramilapath '/.git/refs/heads/master']
fileID = fopen(githashfile,'r');
githash=fgetl(fileID);
fclose(fileID);



% Output settings
setenv('FSLOUTPUTTYPE','NIFTI') % set output type to unarchived .nii

% cfg handling
StdTemplate = cfg.StdTemplate;
TR = cfg.TR;
subj = cfg.inputfolder; disp(['In Data: ' subj]);
subj_out = cfg.outputfolder; disp(['Out Data: ' subj_out]);
if ~exist(subj_out, 'dir')
	mkdir(subj_out);
end
logfilename=[subj_out '/' datestr(now,'yyyymmddHHMMSS') '.log'];
% diary opened although it seems that it does not work with parfor??
diary(logfilename);
diary on
disp(['Log file stored in ' logfilename]); 

% outputs git hash, hardcoded URL
disp(['Pipeline is running bramila version: ' githash]);
disp(['To replicate exact snapshot of pipeline please run:']);
disp(['    git clone git@git.becs.aalto.fi:bml/bramila.git']);
disp(['    git reset --hard ' githash]);
githashoutfile=[subj_out '/githash.txt']; 
fileID=fopen(githashoutfile,'w');
fwrite(fileID,[githash]);
fclose(fileID);



rmvframes = cfg.rmvframes; disp(['How many frames to remove in the beginning: ' num2str(rmvframes)]);
% biopacfile handling
if cfg.drifter == 1
    biopacfile = cfg.biopacfile;
    biopacfile.dtMRI = TR;
    % Check if drifter data is .mat or .acq, assumes there are no other mat
    % files in the folder where the preprocessing will happen
    if ~isempty(dir([subj '/*.acq'])) % if .acq file
        tempfn = dir([subj '/*.acq']);
        biopacfile.name = fullfile(subj,tempfn.name);
        drifterfile=bramila_acq2mat(biopacfile);
        drifterfile=drifterfile{1};
    elseif ~isempty(dir([subj '/*.mat'])) % if .mat file
        tempfn = dir([subj '/*.mat' ]);
        name=ones(1,length(tempfn));
        if length(tempfn)>1
            for i=1:length(tempfn)
                if strcmp(tempfn(i).name(end-6:end), 'cfg.mat') || strcmp(tempfn(i).name(end-6:end),'nii.mat')
                    name(i)=0;
                end
            end
        end
        if sum(name)>1
ut
            disp('detected multiple DRIFTER files (*.mat cfg.mat and *nii.mat are excluded)')
        end
        drifterfile=fullfile(subj,tempfn(name==1).name);
    else
        disp('No DRIFTER files, but drifter==1');
    end
end
% Setup the names
SPGR=sprintf('%s/bet.nii',subj);
SPGRMNI=sprintf('%s/betMNI.nii',subj_out);
inprefix=sprintf('%s/epi',subj);
prefix=sprintf('%s/epi',subj_out);

inFile=sprintf('%s.nii',inprefix);
outSlicedFile=sprintf('%s_st.nii',prefix);
outFileCut=sprintf('%s_cut.nii',prefix);
outFileMCF=sprintf('%s_MCF.nii',prefix);
outFileMean=sprintf('%s_mean.nii',prefix);
outFileBetSingle=sprintf('%s_brain_single.nii',prefix);
outFileBetMask=sprintf('%s_brain_single_mask.nii',prefix);
outFileBet=sprintf('%s_brain.nii',prefix);
outFileMNI=sprintf('%s_STD.nii',prefix);
% outFileSmooth=sprintf('%s_smooth.nii',prefix);

regMatrixStruc2std=sprintf('%s_SPRG2MNI.mat',prefix);
regMatrixFunc2Struc=sprintf('%s_EPI2SPRG.mat',prefix);
regMatrixFunc2std=sprintf('%s_EPI2MNI.mat',prefix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now the main part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Slice timing correction
if isfield(cfg,'sliceseq')
    csvwrite([subj_out '/slices' ],cfg.sliceseq');
    cmd = sprintf('slicetimer -i %s -r %i -v --ocustom=%s/slices -o %s',inFile,TR,subj_out,outSlicedFile);    
else
    %Read in the file header and set the number of slices if slice sequence
    %is not given
	% NOTE - this is only valid for the Siemens scanner since it checks for odd/even slice num.
    tmp_hdr=load_nii_hdr(inFile);
    cfg.slicenum=tmp_hdr.dime.dim(4);
    if mod(cfg.slicenum,2) == 1 %if odd, use default parameters for odd slice number
        cmd = sprintf('slicetimer -i %s --odd -r %i -v -o %s',inFile,TR,outSlicedFile); % -r for TR
    else
        % make the slice order file
        slicelist = [];
        slicelist(1:cfg.slicenum/2) = 2:2:cfg.slicenum;
        slicelist = cat(2,slicelist,1:2:cfg.slicenum);
        csvwrite([subj_out '/slices' ],slicelist');
        cmd = sprintf('slicetimer -i %s -r %i -v --ocustom=%s/slices -o %s',inFile,TR,subj_out,outSlicedFile); % slices is just text file with slice acquisition order
    end
end
system(cmd);
%% Convert to float
cmd = sprintf('fslmaths %s %s -odt float',outSlicedFile,outSlicedFile);
system(cmd);
%% Remove first frames
cmd=sprintf('fslval %s dim4',outSlicedFile); % find out the length of file
[~,tmp2]=system(cmd); vols=str2num(tmp2);
numout=vols - rmvframes;
disp('trying fslroi here');
if exist(outFileCut,'file') == 0  || cfg.overwrite == 1
    fprintf('Removing first %i frames from %s',rmvframes,subj);
    cmd = sprintf('fslroi %s %s %i %i',outSlicedFile,outFileCut,rmvframes,numout);
    disp(cmd);
	system(cmd);
else
    fprintf('Not removing frames from subject %s',subj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: MCFLIRT Motion correction
disp('trying mcflirt')
if exist(outFileMCF,'file') == 0 || cfg.overwrite == 1
    disp('Motion correction');
	cmd=sprintf('mcflirt -in %s -stats -mats -plots -report -rmsrel -rmsabs -o %s',outFileCut,outFileMCF);
    system(cmd);
else
	disp('Skipping motion correction');
end
%Plot MCFLIRT diagnostics
system(sprintf('fsl_tsplot -i %s_MCF.nii.par -t "MCFLIRT estimated rotations (radians)" -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o %s_rot.png',prefix,prefix)); 
system(sprintf('fsl_tsplot -i %s_MCF.nii.par -t "MCFLIRT estimated translations (mm)" -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o %s_trans.png',prefix,prefix));
system(sprintf('fsl_tsplot -i %s_MCF.nii_abs.rms,%s_MCF.nii_rel.rms -t "MCFLIRT estimated mean displacement (mm)" -u 1 -w 640 -h 144 -a absolute,relative -o %s_disp.png',prefix,prefix,prefix)); 			
%BET Brain extraction for the motion corrected functional data
disp('trying bet');
if exist(outFileBetSingle,'file')==0 || cfg.overwrite == 1
    disp('Getting brain mask');
    %Get mean volume
    system(sprintf('fslmaths %s -Tmean %s',outFileMCF,outFileMean));
    system(sprintf('bet %s %s  -f 0.2 -g 0 -o -m -t',outFileMean,outFileBetSingle));
else
    disp('Mask already exists');
end
if exist(outFileBet,'file')==0 || cfg.overwrite == 1
    disp('Masking functional data');
    system(sprintf('fslmaths %s -mas %s %s',outFileMCF,outFileBetMask,outFileBet));
else
    disp('Masked data already exists');
end
%% Step 2.5: Drifter
% Drifter toolbox support to the pipeline was added by Lenka: lenka.vondrackova@aalto.fi
% get biopac measurements filename
% load reference data, for now .acq only?:
if cfg.drifter == 1    
    refdata=load(drifterfile);
    refdata=refdata.refdata;    
    % filtration/normalisation
    if biopacfile.filter == 1
        refdata{1}.filter = 1 ;
        refdata{2}.filter = 1 ;
    else
        refdata{1}.filter = 0 ;
        refdata{2}.filter = 0 ;
    end    
    % cut out the beginning corresponding to rmvtrials
    samptoremove = TR*(refdata{1}.dt^-1)*rmvframes;    
    refdata{1}.data = refdata{1}.data(samptoremove+1:end,1);
    refdata{2}.data = refdata{2}.data(samptoremove+1:end,1);
    refdata{1}.downdt=biopacfile.breath; 
    refdata{2}.downdt=biopacfile.HR;
    refdata{1}.freqlist=biopacfile.freqBreath;
    refdata{2}.freqlist=biopacfile.freqHR;	
    refdata{1}.name=[drifterfile(1:end-4),'_1'];
    refdata{2}.name=[drifterfile(1:end-4),'_2'];
    refdata{1}.outfolder=subj_out;
    refdata{2}.outfolder=subj_out;
    fn=outFileBet;
    nii=load_untouch_nii(fn);
    data.data=double(nii.img);
    data.dt=biopacfile.dtMRI;
    data.mode=cfg.driftermode;
    [data,refdata]=drifter(data, refdata);
    nii.img=data.estimate;
    outDriftFile=[fn(1:end-4) '_cleaned'];
    save_untouch_nii(nii, outDriftFile)
	cfg.drifter_refdata=refdata;
    disp(['saved .... ', outDriftFile, '.nii'])
else
    outDriftFile = outFileBet; % If drifter is not running, take sliced file into the next step
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3: Registration to standard space
disp('trying registration');
if exist(outFileMNI,'file') == 0  || cfg.overwrite == 1
    disp('Registering to MNI standard space');
    % T1 to Standard
    system(sprintf('flirt -in %s -ref %s -omat %s -bins 256 -cost corratio -searchrx -120 120 -searchry -120 120 -searchrz -120 120 -dof 12',SPGR,StdTemplate,regMatrixStruc2std));
    % Func to T1
    system(sprintf('flirt -in %s -ref %s -omat %s -bins 256 -cost corratio -searchrx -120 120 -searchry -120 120 -searchrz -120 120 -dof 9 ',outFileBetSingle,SPGR,regMatrixFunc2Struc));
    % Concatenate to get Func to Standard
    system(sprintf('convert_xfm -concat %s -omat %s %s',regMatrixStruc2std,regMatrixFunc2std,regMatrixFunc2Struc));
    % Transform
    system(sprintf('flirt -in %s -applyxfm -init %s -out %s -paddingsize 0.0 -interp trilinear -ref %s',outDriftFile,regMatrixFunc2std,outFileMNI,StdTemplate));
    system(sprintf('flirt -in %s -applyxfm -init %s -out %s -paddingsize 0.0 -interp trilinear -ref %s',SPGR,regMatrixStruc2std,SPGRMNI,StdTemplate));
else
    disp('Registered data already exists and will not be overwritten');
end
%% Step 4: Bramila
cfg.infile = outFileMNI;
cfg.motionparam = sprintf('%s_MCF.nii.par',prefix);
cfg.save_path = sprintf('%s',subj_out);
[cfg,~] = bramila_clean_signal(cfg);
[cfg,~] = bramila_diagnostics(cfg);
% after the cleaning let's make sure that we use the value in cfg.vol and not cfg.infile
rmfield(cfg,'infile');

%% Step 5: Temporal filtering
if cfg.do_temporal_filtering == 0
	disp('skipping temporal filtering')
else
	cfg=bramila_filter(cfg);
end

%% Step 6: Spatial filtering
if cfg.do_spatial_smooth == 0
	disp('Skipping spatial smoothing')
else
	disp('trying smoothing');
	cfg = bramila_smooth(cfg);
end
%% Write to epi_preprocessed.nii
if ~isfield(cfg,'outfile') || isempty(cfg.outfile)
	cfg.outfile = bramila_savevolume(cfg,cfg.vol,'temprorary file for finishing preprocessing','EPI_tempfile.nii');
else
	disp(['Final preprocessed data from file '  cfg.outfile])
end
curroutfile=cfg.outfile;
target = load_nii(curroutfile);
save_nii(target,[prefix '_preprocessed.nii']);
% we don't need to store the volume
cfg.vol=[];
save(sprintf('%s/bramila/cfg.mat',subj_out),'cfg','-v7.3');
%exit; % importat to exit non-interactive matlab to release the host in the cluster
