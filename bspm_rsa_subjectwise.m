function [rho, fishz] = bspm_rsa_subjectwise(maps, mask, eigenplot, rdmflag)
% 
% 
% USAGE: bspm_rsa_subjectwise(maps, mask, eigenplot, rdmflag)
%
%  ARGUMENTS
%
%

% ----------------------------- Copyright (C) 2014 -----------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 4, rdmflag = 0; end
if nargin < 3, eigenplot = 1; end
if nargin < 2, mask = bspm_greymask; end
if nargin < 1, error('USAGE: bspm_rsa_subjectwise(maps, mask, eigenplot, rdmflag)'); end

% %% MASK
% if iscell(mask), mask = char(mask); end
% mask = bob_reslice(mask,maps{1},1,1);

%% MAPS
if ischar(maps), maps = cellstr(maps); end
all = bspm_read_vol(char(maps), 'reshape', 'implicitmask', 'mask', mask);
% nmaps = length(maps);
% nvox = sum(mask(:) > 0);
% all = zeros(nvox, nmaps);
% for i = 1:nmaps
%     d = bspm_read_vol(maps{i});
%     d = d(:);
%     all(:,i) = d(mask(:) > 0);
% end
all(nanmean(all')==0,:) = 0;
all(find(sum(all'==0)),:) = [];
all(nanmean(isnan(all),2)>0,:) = [];

%% MDS
rho = corr(all, 'rows', 'pairwise');
D = 1 - rho;
[Y,eigvals] = cmdscale(D);
if eigenplot
    figure('color','white');
    plot(1:length(eigvals),eigvals,'bo-');
    if feature('HGUsingMATLABClasses')
        cl = specgraphhelper('createConstantLineUsingMATLABClasses','LineStyle',...
            ':','Color',[.7 .7 .7],'Parent',gca);
        cl.Value = 0;
    else
        graph2d.constantline(0,'LineStyle',':','Color',[.7 .7 .7]);
    end
    axis([1,length(eigvals),min(eigvals),max(eigvals)*1.1]);
    xlabel('Eigenvalue number');
    ylabel('Eigenvalue');
end

%% REPRESENTATIONAL DISSIMILARITY MATRIX OR CORR MAP
figure('color','white');
if rdmflag
    imagesc(D);
else
    imagesc(rho);
end
colormap(jet);
colorbar('SouthOutside');
axis('square');
box off;

%% CLEANUP DIAGONAL
fishz = rho;
fishz(:) = fisherz(rho);
fishz(isinf(fishz)) = NaN;




 
 
 
 
