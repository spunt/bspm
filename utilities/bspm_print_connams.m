function connams = bspm_print_connams(spmmat)
% BSPM_PRINT_CONNAMES
%
%   USAGE: connam = bspm_print_connams(spmmat)
%
if nargin < 1
    spmmat = files('SPM*mat');
    if isempty(spmmat), disp('USAGE: connam = bspm_print_connams(spmmat)'); return; end
end
if iscell(spmmat), spmmat = char(spmmat); end
load(spmmat);
cellprint({SPM.xCon.name}');
if nargin>0, connams = {SPM.xCon.name}'; end 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
