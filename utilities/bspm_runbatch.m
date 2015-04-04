function bspm_runbatch(job)
% USAGE: bspm_runbatch(job)
%
if nargin<1, disp('USAGE: bspm_runbatch(job)'); return; end
spm_jobman('initcfg');
numerrs = 0; 
for i = 1:length(job)
    printmsg(sprintf('Working on Job %d of %d', i, length(job))); 
    try
        spm_jobman('run',job(i));
    catch lasterr
        printmsg(lasterr.message, 'ERROR', '!'); 
        numerrs = numerrs + 1;
        errdata(numerrs).job        = job(i); 
        errdata(numerrs).identifier = lasterr.identifier; 
        errdata(numerrs).message    = lasterr.message; 
        errdata(numerrs).stack      = lasterro.stack;
        continue
    end
end
if numerrs==0
    printmsg('No Errors. Great Job!', 'msgtitle', sprintf('Completed %s', bspm_timestamp), 'msgtop', '='); 
else
    outfile = sprintf('ErrorData_bspm_runbatch_%s.mat', bspm_timestamp); 
    save(outfile, 'errdata'); 
    printmsg(sprintf('Error Data Saved To: %s', outfile), 'msgtitle', sprintf('%d Failed Jobs', numerrs), 'msgtop', '='); 
end
end

 
 
 
 
