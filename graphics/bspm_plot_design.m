function bspm_plot_design(spmmat, rptag, trantag, xlabels)
% BSPM_PLOT_DESIGN
%
%   USAGE: bspm_plot_design(spmmat, rptag, trantag, xlabels)
%
%   INPUTS:
%       spmmat = SPM.mat containing the design
%       rptag = tag to include motion regressors (default = 1)
%       trantag = tag to transpose (default = 0)
%
% Written by Bob Spunt, February 17, 2013

% --------------------- Copyright (C) 2014 ---------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<4, xlabels = []; end
if nargin<3, trantag = 0; end
if nargin<2, rptag = 1; end
if nargin<1, disp('USAGE: bspm_plot_design(spmmat, rptag, trantag, xlabels)'); return; end
if iscell(spmmat), spmmat = char(spmmat); end

% get X matrix and names
load(spmmat)
xmatrix = SPM.xX.xKXs.X; % filtered and whitened design matrix
xname = SPM.xX.name; % regressor names
rowidx = [];
for i = 1:length(xname);
    tmp = xname{i};
    runidx(i) = str2num(tmp(4));
    for r = 1:length(SPM.xX.K)
        rowidx{r} = SPM.xX.K(r).row;
    end
    tmp = tmp(7:end);
    condidx(i) = ~isempty(strfind(tmp,'*bf'));
    if condidx(i)
        tmp = tmp(1:end-6);
    end
    xname{i} = tmp;
end

if ~isempty(xlabels),
    if ischar(xlabels), xlabels = cellstr(xlabels); end
    xname(1:length(xlabels)) = xlabels; 
end

% get rid of unwanted columns
if ~rptag
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
    
    cx = xmatrix(rowidx{r},runidx==r);
    cn = strtrim(xname(runidx==r));
    cn = regexprep(cn,'_',' ');

    % plot
    subplot(1,length(runs),r);
    
    if trantag
          
        imagesc(bspm_scaledata(cx)'); colormap('gray');
        set(gca, 'FontName', 'Arial');
        set(gca, 'YTick', 1:length(cn), 'XTick', (0:25:size(cx,1)));
        set(gca, 'YTickLabel', cn);

    else
        
        imagesc(bspm_scaledata(cx)); colormap('gray');
        
        set(gca, 'FontName', 'Arial');
        set(gca, 'XTick', 1:length(cn), 'YTick', (0:25:size(cx,1)));
        set(gca, 'XTickLabel', '');
        yLim = get(gca,'YLim');
        XTick = get(gca,'XTick');
        y = repmat(yLim(2)+1, length(XTick), 1);
        fontsize = get(gca, 'FontSize');
        hText = text(XTick, y, cn, 'FontSize', floor(.75*(fontsize)));
        set(hText, 'Rotation', 90, 'HorizontalAlignment', 'Right');

    end
    fontsize = get(gca, 'FontSize');
    title(sprintf('Run %d', r), 'FontSize', ceil(1.25*(fontsize)));
    
end




 
 
 
 
 
 
 
 
