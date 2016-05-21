function arrayset(harray, propname, propvalue) 
% ARRAYGET Set property values for array of handles
%
% USAGE: arrayset(harray, propname, propvalue) 
%
% ==============================================
if nargin<2, error('USAGE: arrayset(harray, propname, propvalue) '); end
if size(harray, 1)==1, harray = harray'; end
if ~iscell(propvalue)
    arrayfun(@set, harray, repmat({propname}, length(harray), 1), ...
            repmat({propvalue}, length(harray), 1)); 
else
    if size(propvalue, 1)==1, propvalue = propvalue'; end
    arrayfun(@set, harray, repmat({propname}, length(harray), 1), propvalue); 
end
end