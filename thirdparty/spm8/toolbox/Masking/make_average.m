function Va = make_average(files, out, expr)
% make_average(files, out) -- generate mean image
%
% Email ged.ridgway@gmail.com with any questions

[spm5 select] = check_spm_setup;

% Check if called from batch system
if nargin == 1 && isstruct(files)
    job = files;
    files = job.innames;
    [pth nam ext num] = spm_fileparts(job.outname);
    if isempty(pth)
        % {''} seems to end up as '' in SPM8, hence left hand side of ||
        if ~iscell(job.outdir) || isempty(job.outdir{1})
            pth = pwd;
        else
            pth = job.outdir{1};
        end
    end
    out = fullfile(pth, [nam ext num]);
    expr = job.avgexpr;
end

% Get filenames or SPM.mat
if ~exist('files', 'var') || isempty(files)
    files = select('Select images/masks or SPM.mat'); drawnow
end
if iscellstr(files)
    files = char(files);
end
if size(files, 1) == 1 && ~isempty(regexp(files, 'SPM.mat$', 'once'))
    load(files); % load SPM.mat
    files = char(SPM.xY.P);
end
if size(files, 1) < 1, error('No files specified'); end
Vi = spm_vol(files);

if ~exist('out', 'var') || isempty(out)
    out = 'avg'; % (create in CWD)
end
if isempty(regexp(out, '\.(nii|img)$', 'once'))
    out = [out '.img'];
end

if ~exist('expr', 'var') || isempty(expr)
    expr = 'mean(X)';
end

interp = 1; % linear, though typically identical dimensions anyway
flags = {true, false, interp}; % (dmtx, mask, hold)

Vo = struct(...
    'fname',   out,         ...
    'mat',     Vi(1).mat,   ...
    'descrip', 'mean image' ...
    );
if spm5
    Vo.dim = Vi(1).dim(1:3);
    Vo.dt = [spm_type('float32') Vi(1).dt(2)];
else
    Vo.dim = [Vi(1).dim(1:3) spm_type('float')];
end

Va = spm_imcalc(Vi, Vo, expr, flags);

if exist('job', 'var') 
    % make output compatible with SPM8 dependency, based on spm_cfg_imcalc
    % and its vout subfunction. See also make_average_vout
    clear Va;
    Va.files{1} = out;
end
