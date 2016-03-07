function bnii_4dto3d(fname4d, varargin)
% BNII_4DTO3D
%
%  USAGE: bnii_3dto4d(fname3d, varargin)
% __________________________________________________________________________
%

% ---------------------- Copyright (C) 2015 Bob Spunt ----------------------
%	Created:  2015-04-03
%	Email:     spunt@caltech.edu
% __________________________________________________________________________
def = { ... 
    'basefname',        [], ...
    'delete4d',         0,  ...
	};
vals = setargs(def, varargin);
if nargin==0, mfile_showhelp; fprintf('\t| - VARARGIN DEFAULTS - |\n'); disp(vals); return; end
if iscell(fname4d), fname4d = char(fname4d); end
[basepath, n] = fileparts(fname4d); 
if isempty(basefname)
    basefname = regexprep(n, '\w4D', '', 'ignorecase');  
end
nii = nii_tool('load', fname4d); 
fprintf('| Saving %d 3D images from: %s', size(nii.img, 4), fname4d);
nii_tool('save', nii, fullfile(basepath, strcat(basefname, '.nii')), 1); 
if exist(fname4d, 'file') && delete4d
    fprintf('\n| Deleting 4D image'); 
    rm(fname4d, 'verbose', 0); 
end
fprintf('\n - Process Complete -\n');
end
% - SUBFUNCTIONS
function argstruct = setargs(defaultargs, varargs)
    if nargin < 1, mfile_showhelp; return; end
    if nargin < 2, varargs = []; end
    defaultargs = reshape(defaultargs, 2, length(defaultargs)/2)'; 
    if ~isempty(varargs)
        if mod(length(varargs), 2)
            error('Optional inputs must be entered as "''Name'', Value" pairs, e.g., myfunction(''arg1'', val1, ''arg2'', val2)'); 
        end
        arg = reshape(varargs, 2, length(varargs)/2)';
        for i = 1:size(arg,1)
           idx = strncmpi(defaultargs(:,1), arg{i,1}, length(arg{i,1}));
           if sum(idx) > 1
               error(['Input "%s" matches multiple valid inputs:' repmat('  %s', 1, sum(idx))], arg{i,1}, defaultargs{idx, 1});
           elseif ~any(idx)
               error('Input "%s" does not match a valid input.', arg{i,1});
           else
               defaultargs{idx,2} = arg{i,2};
           end
        end
    end
    for i = 1:size(defaultargs,1), assignin('caller', defaultargs{i,1}, defaultargs{i,2}); end
    if nargout>0, argstruct = cell2struct(defaultargs(:,2), defaultargs(:,1)); end
end
function mfile_showhelp(varargin)
    ST = dbstack('-completenames');
    if isempty(ST), fprintf('\nYou must call this within a function\n\n'); return; end
    eval(sprintf('help %s', ST(2).file));  
end