function [xmatrix, xname] = bspm_check_design(spmmat, r_tag)
% BSPM_CHECK_DESIGN
%
%   USAGE: [xmatrix xname] = bspm_check_design(spmmat, r_tag)
%
%   INPUTS:
%       spmmat = SPM.mat containing the design
%       r_tag = tag to include motion regressors
%
%   OUTPUTS:
%       xmatrix = filtered and whitened design matrix
%       xname = regressor names
%
% Written by Bob Spunt, February 17, 2013

% --------------------- Copyright (C) 2014 ---------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<2, mfile_showhelp; return; end
if iscell(spmmat), spmmat = char(spmmat); end

% get X matrix and names
load(spmmat)
xmatrix = SPM.xX.xKXs.X; % filtered and whitened design matrix
xname = SPM.xX.name; % regressor names
for i = 1:length(xname);
    tmp = xname{i};
    runidx(i) = str2num(tmp(4));
    tmp = tmp(7:end);
    condidx(i) = ~isempty(strfind(tmp,'*bf'));
    if condidx(i)
        tmp = tmp(1:end-6);
    end
    xname{i} = tmp;
end

% get rid of unwanted columns
if ~r_tag
    xmatrix = xmatrix(:,condidx);
    xname = xname(condidx);
    runidx = runidx(condidx);
else
    idx = strcmp(xname,'constant');
    xmatrix(:,idx) = [];
    xname(idx) = [];
    runidx(idx) = [];
end

runs = unique(runidx);
figure('Color','white')
for r = runs
    
    cx = xmatrix(:,runidx==r);
    cn = xname(runidx==r);
    rmatrix = corr(cx);
    fprintf('\nRUN %d\n',r);
    disptable(rmatrix, cn, cn, '%.2f');
    
    % plot
    subplot(1,length(runs),r); 
    imagesc(rmatrix); colormap('gray');
    title(sprintf('Run %d', r));
    
%     [x y] = meshgrid(1:length(cn));
%     textStrings = num2str(rmatrix(:), '%0.2f');
%     textStrings = strtrim(cellstr(textStrings));
%     hStrings = text(x(:),y(:),textStrings(:), 'HorizontalAlignment', 'center');
%     midValue = mean(get(gca, 'CLim'));
%     textColors = repmat(rmatrix(:) > midValue,1,3);
%     set(hStrings,{'Color'},num2cell(textColors,2));
    
    set(gca, 'XTick', 1:length(cn), 'YTick', 1:length(cn));
    set(gca, 'XTickLabel', cn, 'YTickLabel', cn);
    
    fprintf('Intercorrelation Summary\n');
    measures = {'MEAN' 'MAX' 'MIN'};
    for c = 1:length(cn)
        tmp = rmatrix(:,c);
        tmp(c) = [];
        summary(c,1) = nanmean(tmp);
    end
    for c = 1:length(cn)
        tmp = rmatrix(:,c);
        tmp(c) = [];
        summary(c,2) = max(tmp);
    end
    for c = 1:length(cn)
        tmp = rmatrix(:,c);
        tmp(c) = [];
        summary(c,3) = min(tmp);
    end
    disptable(summary, measures, cn, '%.2f');
    
    
end





 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
