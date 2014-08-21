function bspm_check_design_rank(spmmat)
% BSPM_CHECK_DESIGN_RANK
%
%   USAGE: bspm_check_design_rank(spmmat)
%
%   INPUTS:
%       spmmat = SPM.mat containing the design
%
% Written by Bob Spunt, February 17, 2013

% ------------ Copyright (C) 2014 ------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, disp('bspm_check_design_rank(spmmat)'); return; end
if iscell(spmmat), spmmat = char(spmmat); end

% get X matrix and names
load(spmmat)
xmatrix = SPM.xX.xKXs.X; % filtered and whitened design matrix
xname = SPM.xX.name; % regressor names
est = spm_SpUtil('IsCon', SPM.xX.xKXs);
if sum(est==0)>0
    
    fprintf('\nPROBLEM\n');
    
    for i = 1:length(xname);
        tmp = xname{i};
        runidx(i) = str2num(tmp(4));
        tmp = tmp(7:end);
        condidx(i) = ~isempty(strfind(tmp,'*bf'));
        if condidx(i)
            tmp = tmp(1:end-6);
        end
        if est(i)==0
            fprintf('Run %d\t%s\n', runidx(i), tmp);
        end
     
    end
    
else
    
    fprintf('\nAll Estimable\n');

end
    





 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
