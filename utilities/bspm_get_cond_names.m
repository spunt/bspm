function [condnames runidx] = bspm_get_cond_names(spmmat)
% --------------------- Copyright (C) 2014 ---------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<1, spmmat = 'SPM.mat'; end
if iscell(spmmat), spmmat = char(spmmat); end

% get X matrix and names
load(spmmat)
condnames = SPM.xX.name'; % regressor names

for i = 1:length(condnames);
    tmp = condnames{i};
    runidx(i) = str2num(tmp(4));
    tmp = tmp(7:end);
    condidx(i) = ~isempty(strfind(tmp,'*bf'));
    if condidx(i)
        tmp = tmp(1:end-6);
    end
    condnames{i} = tmp;
end

% create sessidx
runidx = zeros(size(condnames,1),1);
for r = 1:length(SPM.Sess)
    runidx(SPM.Sess(r).row,1) = r;
end
    




 
 
 
 
