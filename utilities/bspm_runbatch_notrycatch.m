function bspm_runbatch_notrycatch(job, doparfor)
% USAGE: bspm_runbatch(job, doparfor)
%
if nargin<1, mfile_showhelp; return; end
if nargin<2, doparfor = 0; end
% | Initialize SPM config 
bspm_init;
njob = length(job);
if doparfor
    parfor i = 1:njob
        printmsg(sprintf('Working on Job %d of %d', i, length(job))); 
        spm_jobman('run',job(i));
    end
else
    for i = 1:njob
        printmsg(sprintf('Working on Job %d of %d', i, length(job))); 
        spm_jobman('run',job(i));
    end
end
end

 
 
 
 
