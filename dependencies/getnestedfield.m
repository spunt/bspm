function F = getnestedfield(S, EXP)
% GETNESTEDFIELD  Find a field in a structure and return its contents
% 
%  USAGE: F = getnestedfield(S, EXP)
% __________________________________________________________________________
%  INPUTS
% 	S       The struct to search  
% 	EXP     Pattern used to find field (passed to REGEXP)
% __________________________________________________________________________
%  OUTPUT
% 	F       Field contents.
%           If no fields are found, is an empty array. 
%           If multiple fields found, is a N x 2 cell array where N is the
%           number of matching fields, first column contains full locations
%           of each match, and second column contains their values (this is
%           raw output from NSTRUCT2CELL, see CREDIT below).
% __________________________________________________________________________
%  CREDIT
% 	Includes MATLAB File Exchange contribution by Mathias Benedek: 
%   www.mathworks.com/matlabcentral/fileexchange/29908-nstruct2cell
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-03-17
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 2, mfile_showhelp; return; end
F   = [];               % init output
C   = nstruct2cell(S);  % call NSTRUCT2CELL for fields/contents
CF  = cellfun(@(x) x(end), regexp(C(:,1), '\.', 'split')); % get fieldnames
IDX = find(~cellfun('isempty', regexp(CF, EXP))); % find indices of EXP

% | if no field is found matching EXP...
if isempty(IDX), disp('No fields found matching that expression'); return; end

% | if more than one field is found matching EXP...
if length(IDX) > 1, fprintf('\nMultiple fields matching that value:\n'); disp(C(IDX,:)); F = C(IDX,:); return; end
    
% | if the temperature of the porridge is just right...
F = C{IDX,2};

end
% | SUBFUNCTIONS
function C = nstruct2cell(S, BRANCH)
%C = NSTRUCT2CELL(S)
% Recursive function that converts a nested struct S with a total of n sub-fields into a nx2 cell array C.
% The first column of C lists the full names of the sub-fields and the second column contains the respective content.
% Mathias Benedek 201-01-04
if nargin == 1  %First level of struct
    BRANCH = inputname(1);
end
C = {};
if ~isstruct(S)     % End of struct-branch, no further fields: read content
    C = {BRANCH, S};
else
    n = numel(S);
    if n == 1   % (non-array) struct: parse fields
        fn = fieldnames(S);
        for ii = 1:length(fn)
            C = [C; nstruct2cell(S.(fn{ii}), [BRANCH,'.',fn{ii}])];
        end

    else        % struct-array: parse array elements
        for jj = 1:n
            C = [C; nstruct2cell(S(jj), [BRANCH,'(',num2str(jj),')'])];
        end

    end

end
end
function mfile_showhelp(varargin)
    ST = dbstack('-completenames');
    if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
    eval(sprintf('help %s', ST(2).file));  
end
