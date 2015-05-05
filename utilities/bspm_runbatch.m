function bspm_runbatch(job, doparfor)
% USAGE: bspm_runbatch(job, doparfor)
%
if nargin<1, disp('USAGE: bspm_runbatch(job, doparfor)'); return; end
if nargin<2, doparfor = 0; end
% | Initialize SPM config 
bspm_init;
njob = length(job);
[joblog(1:njob).job]    = deal(job{:});
[joblog(1:njob).status] = deal('Success');
[joblog(1:njob).errdata] = deal([]); 
numerrs = 0; 
if doparfor
    parfor i = 1:njob
        printmsg(sprintf('Working on Job %d of %d', i, length(job))); 
        try
            spm_jobman('run',job(i));
        catch lasterr
            disp(lasterr.message);
            joblog(i).status = 'Error';  
            joblog(i).errdata = lasterr; 
            numerrs = numerrs + 1;
            continue
        end
    end
else
    for i = 1:njob
        printmsg(sprintf('Working on Job %d of %d', i, length(job))); 
        try
            spm_jobman('run',job(i));
        catch lasterr
            disp(lasterr.message); 
            joblog(i).status = 'Error';  
            joblog(i).errdata = lasterr; 
            numerrs = numerrs + 1;
            continue
        end
    end
end
if numerrs==0
    printmsg('No Errors. Great Job!', 'msgtitle', sprintf('Completed %s', bspm_timestamp), 'msgtop', '=', 'msgbottom', '='); 
else
    outfile = sprintf('JOBLOG_BSPM_RUNBATCH_%s.mat', bspm_timestamp); 
    save(outfile, 'joblog'); 
    printmsg(sprintf('Error Data Saved To: %s', outfile), 'msgtitle', sprintf('%d Failed Jobs', numerrs), 'msgtop', '=', 'msgbottom', '='); 
end
end

 
 
 
 
