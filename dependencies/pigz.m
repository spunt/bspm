function pigz(in, varargin)
% PIGZ Use system PIGZ to de/compress
%
%  USAGE: pigz(in, varargin)
% __________________________________________________________________________
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-03
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
        'forceoverwrite',   1,  ...
        'keeporiginal',     0,  ...
        'verbose',          0,  ...
        };
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if ischar(in)
    if regexp(in, '\*')
        in = files(in);
        if isempty(in)
            disp('No files found with that wildcard pattern');
            return; 
        end
    else
        in = cellstr(in); 
    end
end
[fp,fn,fe]  = fileparts(in{1}); 
flgidx      = find([strcmp(fe, '.gz') verbose forceoverwrite keeporiginal]);
if ~isempty(flgidx)
    flgstr  = {'d' 'v' 'f' 'k'};
    flags   = strcat('-', flgstr(flgidx), {' '});
    flags   = strtrim([flags{:}]);
    basecmd     = sprintf('pigz %s', flags);
else
    basecmd     = 'pigz ';
end

fprintf('\nRunning %s on %d files\n', basecmd, length(in)); 
for i = 1:length(in)
    fprintf('%s: %s\n', printcount(i, length(in)), in{i}); 
    [err, str] = system(sprintf('%s "%s"', basecmd, in{i})); 
    if err, fprintf('Error: %s\n', str); end
end
