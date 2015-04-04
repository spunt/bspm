function bspm_runbatch_cluster(job)
% USAGE: bspm_runbatch_cluster(job)
%
if nargin<1, disp('USAGE: bspm_runbatch(job)'); return; end
spm_jobman('initcfg');
for i = 1:length(job), spm_jobman('run', job{i}); end
end

 
 
 
 
