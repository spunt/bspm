function [status, msg] = checkfields(struct, fnames)
% CHECKFIELDS(STRUCT, FIELDNAMES)
%
% -----------------------------------------------------------
if nargin < 2, error('CHECKFIELDS(STRUCT, FIELDNAMES)'); end
isnotidx = ~isfield(struct, fnames);
if any(isnotidx)
    status = 0; 
    msg = sprintf(['\nINPUT STRUCTURE IS MISSING FIELDS:\n' repmat('  %s\n',1,sum(isnotidx))], fnames{isnotidx});
else
    status = 1;
    msg = sprintf('\nInput looks OK\n');
end

