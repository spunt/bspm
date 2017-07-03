function logtransform = tbx_cfg_logtransform
% Configuration file for toolbox 'Logtransform'
%_______________________________________________________________________
% Copyright (C) 2013 Dave Langers <dave@langers.nl>

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','Logtransform')); end

% ---------------------------------------------------------------------
% data Session
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Session';
data.help    = {'Select images for this session.'};
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];
% ---------------------------------------------------------------------
% generic Data
% ---------------------------------------------------------------------
generic        = cfg_repeat;
generic.tag    = 'generic';
generic.name   = 'Data';
generic.help   = {'Add new image session. All images in a session should have similar contrast; different sessions may belong to the same or to different subjects.'};
generic.values = {data };
generic.num    = [1 Inf];
%--------------------------------------------------------------------------
% ascale Automatic scaling
%--------------------------------------------------------------------------
ascale      = cfg_const;
ascale.tag  = 'ascale';
ascale.name = 'Automatic scaling';
ascale.val  = {1};
ascale.help = {'A suitable value of the scaling constant X0 in the transformation Y = 100*ln(X/X0) will be automatically determined on the basis of the first image of each session. The output images will retain a natural tissue contrast.'};
%--------------------------------------------------------------------------
% val Value
%--------------------------------------------------------------------------
val         = cfg_entry;
val.tag     = 'val';
val.name    = 'Value';
val.help    = {'Enter a constant value for the scaling constant X0 in the transformation Y = 100*ln(X/X0); the same value will be applied to all sessions.'};
val.strtype = 'e';
val.num     = [1 1];
val.val     = {1};
%--------------------------------------------------------------------------
% vscale Scale by value
%--------------------------------------------------------------------------
vscale      = cfg_branch;
vscale.tag  = 'vscale';
vscale.name = 'Scale by value';
vscale.val  = {val };
vscale.help = {'Enter a constant value for the scaling constant X0 in the transformation Y = 100*ln(X/X0); the same value will be applied to all sessions.'}';
%--------------------------------------------------------------------------
% img Image
%--------------------------------------------------------------------------
img         = cfg_files;
img.tag     = 'img';
img.name    = 'Image';
img.help    = {'Specify an image that contains voxel-wise values of the scaling constant X0 in the transformation Y = 100*ln(X/X0); the same image will be applied to all sessions. All input images must be in the same coordinate space and have the same dimensions as this image.'};
img.filter  = {'image'};
img.ufilter = '.*';
img.num     = [1 1];
%--------------------------------------------------------------------------
% iscale Scale by image
%--------------------------------------------------------------------------
iscale      = cfg_branch;
iscale.tag  = 'iscale';
iscale.name = 'Scale by image';
iscale.val  = {img };
iscale.help = {'Specify an image that contains voxel-wise values of the scaling constant X0 in the transformation Y = 100*ln(X/X0); the same image will be applied to all sessions. All input images must be in the same coordinate space and have the same dimensions as this image.'}';
%--------------------------------------------------------------------------
% scaling Scaling
%--------------------------------------------------------------------------
scaling        = cfg_choice;
scaling.tag    = 'scaling';
scaling.name   = 'Scaling';
scaling.val    = {ascale };
scaling.help   = {'The value of the scaling constant X0 in the transformation Y = 100*ln(X/X0) can be determined automatically, or be manually specified as a constant value or an image.'};
scaling.values = {ascale vscale iscale };
% ---------------------------------------------------------------------
% clipneg Negative values
% ---------------------------------------------------------------------
clipneg        = cfg_menu;
clipneg.tag    = 'clipneg';
clipneg.name   = 'Negative values';
clipneg.help   = {'Negative values that may occur after the transformation can either be retained, or clipped to zero.'};
clipneg.labels = {
                 'Retain'
                 'Clip to zero'
}';
clipneg.values  = {false true};
clipneg.val     = {false};
% ---------------------------------------------------------------------
% dtype Data type
% ---------------------------------------------------------------------
dtype        = cfg_menu;
dtype.tag    = 'dtype';
dtype.name   = 'Data type';
dtype.help   = {'Data-type of output images. SAME indicates the same datatype as the images from the original session.'};
dtype.labels = {
               'SAME'
               'UINT8   - unsigned 8-bit integer'
               'UINT16  - unsigned 16-bit integer'
               'UINT32  - unsigned 32-bit integer'
               'INT8    - signed 8-bit integer'
               'INT16   - signed 16-bit integer'
               'INT32   - signed 32-bit integer'
               'FLOAT32 - single precision floating point'
               'FLOAT64 - double precision floating point'
}';
dtype.values  = {0 spm_type('uint8') spm_type('uint16') spm_type('uint32') spm_type('int8') spm_type('int16') spm_type('int32') spm_type('float32') spm_type('float64')};
dtype.val     = {0};
%--------------------------------------------------------------------------
% Log Transform
%--------------------------------------------------------------------------
logtransform      = cfg_exbranch;
logtransform.tag  = 'logtransform';
logtransform.name = 'Log-transform';
logtransform.val  = {generic scaling clipneg dtype};
logtransform.help = {
                    'This routine applies a logarithmic transformation Y = 100*ln(X/X0) to series of images. This transforms small relative deviations from baseline into absolute signals. By applying this transformation during preprocessing, the BOLD responses that are extracted using the GLM are automatically expressed as a percentage signal change, as opposed to arbitrary original scanner units.'
                    ''
                    'When applying this transformation, the ''Global normalisation'' option in the fMRI model specification should be assigned default ''None''.'
}';
logtransform.prog = @tbx_run_logtransform;
logtransform.vout = @vout_logtransform;

%==========================================================================
function dep = vout_logtransform(job)
for k=1:numel(job.data)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Log-transformed Images (Session %d)', k);
    dep(k).src_output = substruct('.','session', '()',{k}, '.','lfiles');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;

