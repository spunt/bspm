function bspm_setdefaults(def, args)
% SETDEFAULTS Set default input arguments
%
%  USAGE: USAGE: setdefaults(def, args)  
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-03-11
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 1, disp('USAGE: setdefaults(def, args)'); return; end
if nargin < 2, args = []; end
def = reshape(def, 2, length(def)/2)'; 
if ~isempty(args)
    arg = reshape(args, 2, length(args)/2)';
    for i = 1:size(arg,1)
       idx = strncmpi(def(:,1), arg{i,1}, length(arg{i,1})); 
       if any(idx), def{idx, 2} = arg{i, 2}; end
    end
end
for i = 1:size(def,1), assignin('caller', def{i,1}, def{i,2}); end
end