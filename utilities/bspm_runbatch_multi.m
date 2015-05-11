function bspm_runbatch_multi(job)
% USAGE: bspm_runbatch(job)
%
if nargin<1, disp('USAGE: bspm_runbatch(job)'); return; end
% | Initialize SPM config 
bspm_init;
njob = length(job);
for i = 1:njob
    try
        spm_jobman('run',job(i));
    catch lasterr
        rethrow(lasterr); 
    end
end

end

 
 
 
 
