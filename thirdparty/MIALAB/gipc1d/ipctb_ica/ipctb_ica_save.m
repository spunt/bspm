function ipctb_ica_save(varargin)
% Function to save the files. Automatically uses the version 6
% compatibility on higher Matlab versions

whichVersion = str2num(version('-release'));

if isempty(whichVersion)
    whichVersion = 14;
end

% loop over nargs
for ii = 2:nargin
    if ~strcmp(varargin{ii}(1), '-')
        % get only the variable names
        eval([varargin{ii}, ' = ', 'evalin(''caller'', varargin{ii});']);
    end
end

% Make the MAT file version 6 compatible
if whichVersion > 13
    save(varargin{1}, varargin{2:end}, '-v6');
else
    save(varargin{1}, varargin{2:end});
end