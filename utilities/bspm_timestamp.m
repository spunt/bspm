function str = bspm_timestamp(dateonly)
% USAGE: str = bspm_timestamp(dateonly)
if nargin==0, dateonly = 0; end
if dateonly
    str = datestr(now,'mm-DD-YYYY'); 
else
    str = sprintf('%s_%s', datestr(now,'mm-DD-YYYY'), lower(strtrim(datestr(now, 'HHMMpm')))); 
end



