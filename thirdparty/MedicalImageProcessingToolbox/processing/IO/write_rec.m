function write_rec(fname,data,format)

% save REC file
% fname  : name of the .rec to be saved
% data   : data to be saved in .rec file
% format : optional, 'uint16' by default
% Claudia Prieto

if (nargin == 2)
  format = 'uint16';
end

pid = fopen(fname,'w');
if (pid == -1),
    error('could not create file...');
end
fwrite(pid,[data(:)],format);   %[data(:); data(:)]
fclose(pid);
