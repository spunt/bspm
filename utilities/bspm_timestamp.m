function str = bspm_timestamp(dateonly, inclsecs)
% BSPM_TIMESTAMP
%   
% USAGE: str = bspm_timestamp(dateonly)
%
if nargin==0, dateonly = 0; end
if nargin<2, inclsecs = 0; end
dayopt      = {datestr(now,'YYYYmmDD') strtrim(datestr(now,'YYYYmmmDD'))};
timeopt     = {lower(strtrim(datestr(now, 'HHMMpm'))) strtrim(datestr(now,'HHMMSSpm'))};
if dateonly
	str = upper(dayopt{1});
else
    str = sprintf('%s_%s', dayopt{1}, timeopt{1 + inclsecs});
end



