function out = tbx_run_logtransform(job)
% Job execution function for toolbox 'Logtransform'
% takes a harvested job data structure and performs a log-transformation.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2013 Dave Langers <dave@langers.nl>

% Dave Langers

spm_progress_bar('Init', sum(cellfun(@numel, job.data)), 'Log-transforming Images');
n = 0;

if isfield(job.scaling, 'vscale')
  scale = job.scaling.vscale.val;
elseif isfield(job.scaling, 'iscale')
  scale = spm_read_vols(spm_vol(job.scaling.iscale.img{1}));
end

for i=1:numel(job.data)
  
  pinfo = zeros(3, 1);
  
  for k = 1:numel(job.data{i})
    
    n = n+1;
    spm_progress_bar('Set',n);
    
    V = spm_vol(job.data{i}{k});
    X = spm_read_vols(V);

    if k == 1 && isfield(job.scaling, 'ascale'), scale = max(mean(X(X > mean(X(:))/10))/10, 0); end

    V.fname = spm_file(V.fname,'prefix','l');
    if job.dtype ~= 0, V.dt(1) = job.dtype; end
    V.pinfo = pinfo;
    
    if job.clipneg
      X(X < scale) = scale;
    else
      X(X <= 0) = nan;
    end
    X = 100*log(X./scale);

    V = spm_write_vol(V, X);
    pinfo = V.pinfo;
    
    out.session(i).lfiles{k} = V.fname;
    
  end
  
end

spm_progress_bar('Clear')
