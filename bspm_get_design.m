function [xmatrix xname sessidx] = bspm_get_design(spmmat, r_tag)
% BSPM_GET_DESIGN  Extract filtered/whitened design matrix from SPM.mat
%
%   USAGE: [xmatrix xname] = bspm_get_design(spmmat, r_tag)
%
%   INPUTS:
%       spmmat = SPM.mat containing the design
%       r_tag = tag to include motion regressors (default=1)
%
%   OUTPUTS:
%       xmatrix = filtered and whitened design matrix
%       xname = regressor names
%       sessidx = session # of rows of xmatrix
%
% Written by Bob Spunt, February 17, 2013

% -------------------------- Copyright (C) 2014 --------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, disp('USAGE: [xmatrix xname sessidx] = bspm_get_design(spmmat, r_tag)'); return; end
if nargin<2, r_tag = 1; end
if iscell(spmmat), spmmat = char(spmmat); end

% get X matrix and names
load(spmmat)
xmatrix = SPM.xX.xKXs.X; % filtered and whitened design matrix
xname = SPM.xX.name; % regressor names

% check for rWLS
flag_rwls = 0;
if isfield(SPM,'ResStats')
    flag_rwls = 1;
end

% get condition names
for i = 1:length(xname);
    tmp = xname{i};
    runidx(i) = str2num(tmp(4));
    tmp = tmp(7:end);
    condidx(i) = ~isempty(strfind(tmp,'*bf'));
    if condidx(i)
        tmp = tmp(1:end-6);
    end
    tmp = regexprep(tmp,'\^1','');
    xname{i} = tmp;
end

% get rid of unwanted columns
if ~r_tag
    xmatrix = xmatrix(:,condidx);
    xname = xname(condidx);
    runidx = runidx(condidx);
else
    if flag_rwls
        xname = regexprep(xname,'constant','rWLS');
    else
        xname( strcmp(xname,'constant')) = [];
    end
end

% create sessidx
sessidx = zeros(size(xmatrix,1),1);
for r = 1:length(SPM.Sess)
    sessidx(SPM.Sess(r).row) = r;
end


 
 
 
 
