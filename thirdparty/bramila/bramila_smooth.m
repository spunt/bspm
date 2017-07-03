function cfg = bramila_smooth(varargin)
if nargin==1
    cfg = varargin{1};
else
    cfg = varargin{1};
    group_cfg = varargin{2};
end
% wrapper function for spatial smoothing

fprintf('Running spatial smoothing\n');
smooth_cfg = [];
smooth_cfg.programpath=cfg.bramilapath;
smooth_cfg.separator=cfg.separator;
smooth_cfg.FSLDIR = cfg.FSLDIR;

delete_EPI_temp = 0;
if nargin > 1 && iscell(cfg) % handle the group case
    if group_cfg.write == 1 || (~isfield(cfg,'vol') || (isfield(cfg,'vol') && isempty(cfg.vol)))
        smooth_cfg.infile = cfg.infile;
    else
        fprintf('..creating tempfile\n');
        smooth_cfg.infile = bramila_savevolume(cfg,cfg.vol,'temprorary file for smoothing','EPI_tempfile.nii');
        delete_EPI_temp = 1;
    end
else
    if (~isfield(cfg,'vol') || (isfield(cfg,'vol') && isempty(cfg.vol)))
        smooth_cfg.infile = cfg.infile;
    else
        fprintf('..creating tempfile\n');
        smooth_cfg.infile = bramila_savevolume(cfg,cfg.vol,'temprorary file for smoothing','EPI_tempfile.nii');
        delete_EPI_temp = 1;
    end
    
end
smooth_cfg.FWHM=cfg.smooth_FWHM;

try
    if strcmp(cfg.smooth_method,'SPM')
        resfile = run_SPM_smooth(smooth_cfg);
    elseif strcmp(cfg.smooth_method,'AFNI')
        resfile = run_AFNI_smooth(smooth_cfg);
    elseif strcmp(cfg.smooth_method,'FSL')
        resfile = run_FSL_smooth(smooth_cfg);
    elseif strcmp(cfg.smooth_method,'FSLgauss')
        resfile = run_FSLgauss_smooth(smooth_cfg);
    else
        error('Unknown smoothing method!')
    end
catch err
    warning('!! Smoothing failed! Skipping... !!');
    rethrow(err);
    return;
end

fprintf('..finalizing files\n');

nii = load_nii(resfile);
vol = nii.img;
cfg.vol=vol;
% Filename is hardcoded now, maybe for future it would make sense to
% make it appended depending on which steps were done
disp(cfg.infile_orig)
disp(resfile)
disp(num2str(nargout))
if cfg.write==1 || nargout<1
    cfg.infile = bramila_savevolume(cfg,vol,'EPI volume after masking, detrending, motion(+tissue) regression, filtering (check if you did) and smoothing','mask_detrend_fullreg_filtered_smoothed.nii');
	cfg.outfile = cfg.infile;

	if ~strcmp(cfg.infile_orig,resfile)
		delete(resfile);
	else
		warning('!! For some reason, temporary smoothed EPI file is the original EPI file (did you skip some steps or mix names?), cannot delete !!')
	end
else
    cfg.outfile=[]
	if(isfield(cfg,'outfile'))
	   cfg = rmfield(cfg,'outfile');
	end
end

if delete_EPI_temp==1
    if ~strcmp(cfg.infile_orig,smooth_cfg.infile)
        delete(smooth_cfg.infile);
    else
        warning('!! For some reason, temporary EPI file is the original EPI file (did you skip some steps?), cannot delete !!')
    end
end

end

function resfile = run_FSL_smooth(cfg)
% Using SUSAN
% fsl_root_command = 'fsl5.0-'; % Triton uses different version of FSL by
% default. One solution could be to have FSLDIR defined in
% check_parameters.m
setenv('FSLOUTPUTTYPE','NIFTI') % set output type to unarchived .nii
infile = cfg.infile;
fslpath=cfg.FSLDIR;
FWHM=cfg.FWHM/2.35; % FSL takes sigma

[orig_path,orig_id,~]=fileparts(infile);
resfile = [orig_path,cfg.separator,orig_id,'_smooth.nii'];

fprintf('..running FSL SUSAN\n');

[a,b]=unix(sprintf('%ssusan %s -1 %f 6 1 0 %s',fslpath,infile,FWHM,resfile),'-echo');
if a~=0
    disp(b);
    error('FSL susan smoothing failed')
end
if exist(resfile,'file')==0
    % nii file not found, we must convert
    str = [fslpath,'fslchfiletype NIFTI ',resfile];
    [a,b] = unix(str);
    if a~=0
        disp(b);
        error('FSL fslchfiletype failed')
    end
end

end

function resfile = run_FSLgauss_smooth(cfg)
setenv('FSLOUTPUTTYPE','NIFTI')
infile = cfg.infile;
FWHM=cfg.FWHM/2.35; % FSL takes sigma
fslpath=cfg.FSLDIR;

[orig_path,orig_id,~]=fileparts(infile);
resfile = [orig_path,cfg.separator,orig_id,'_smooth.nii'];

fprintf('..running FSL maths\n');

[a,b]=unix(sprintf('%sfslmaths %s -kernel gauss %f -fmean %s',fslpath,infile,FWHM,resfile),'-echo');
if a~=0
    disp(b);
    error('FSL maths smoothing failed')
end
if exist(resfile,'file')==0
    % nii file not found, we must convert
    str = [fslpath,'fslchfiletype NIFTI ',resfile];
    [a,b] = unix(str);
    if a~=0
        disp(b);
        error('FSL fslchfiletype failed')
    end
end
end

function output = run_SPM_smooth(cfg)

infile = cfg.infile;

FWHM=cfg.FWHM;

[ total_scan ] = get_nii_frame(infile);

[orig_path,orig_id,ext]=fileparts(infile);
output = [orig_path,filesep,'smoothed_',orig_id,ext];

for i=1:total_scan
    matlabbatch{1}.spm.spatial.smooth.data{i}= [infile,',',num2str(i)];
end
matlabbatch{1}.spm.spatial.smooth.fwhm = FWHM*ones(1,3);
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 'smoothed_';

spm('defaults','fmri');
spm_jobman('initcfg');

fprintf('..running SPM smooth\n');
spm_jobman('run',matlabbatch);

if ~exist(output,'file')
    error('Resultfile not found!')
end

end


function resfile = run_AFNI_smooth(cfg)

infile = cfg.infile;
FWHM=cfg.FWHM;

[orig_path,orig_id,~]=fileparts(infile);
id = [orig_path,cfg.separator,orig_id,'_tempfile'];

command = [' -FWHM ',num2str(FWHM),' -input ',infile,' -prefix ',id];

fprintf('..running AFNI smooth with FWHM=%3.1f. Input file: %s\n',FWHM,infile);

[a,b] = unix(['chmod +x ',cfg.programpath,cfg.separator,'external',cfg.separator,'bin',cfg.separator,'3dBlurToFWHM']);
[a,b] = unix([cfg.programpath,cfg.separator,'external',cfg.separator,'bin',cfg.separator,'3dBlurToFWHM ',command],'-echo');

if a~=0
    disp(b);
    error('AFNI 3dBlurToFWHM failed')
end

if exist([id,'+tlrc.BRIK'],'file')
    id = [id,'+tlrc'];
elseif exist([id,'+orig.BRIK'],'file')
    id = [id,'+orig'];
else
    error('No output files found!');
end

command = ['-prefix ',orig_path,cfg.separator,orig_id,'_smooth ',id];

[a,b] = unix(['chmod +x ',cfg.programpath,cfg.separator,'external',cfg.separator,'bin',cfg.separator,'3dAFNItoNIFTI']);
[a,b] = unix([cfg.programpath,cfg.separator,'external',cfg.separator,'bin',cfg.separator,'3dAFNItoNIFTI ',command]);

if a~=0
    disp(b)
    error('AFNI 3dAFNItoNIFTI failed')
end

if exist([id,'.HEAD'],'file')
    delete([id,'.HEAD']);
end
if exist([id,'.BRIK'],'file')
    delete([id,'.BRIK']);
end

resfile = [orig_path,cfg.separator,orig_id,'_smooth.nii'];
if ~exist(resfile,'file')
    error('Resultfile not found!')
end

end

