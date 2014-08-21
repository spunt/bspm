function result = bspm_imcompare(in, labels, sdcut)
% BSPM_IMCOMPARE  Imcalc without writing image
%
%   USAGE: bspm_imcompare(in, labels, sdcut)
%       
%       in  =  array of images
%       labels = string specifying option
%       sdcut = default is 2.5
%
%           PAIRWISE COMPARISONS (2 ONLY)
%               'meanSS'        - mean sum of squared differences
%               'corr'          - correlation across non-zero voxels
%               
%           DISTANCE FROM MEAN IMAGE (2 OR MORE)
%               'MmeanSS'       - mean SS from mean across images
%               'Mcorr'         - correlation with mean image
%
% CREATED: Bob Spunt, Ph.D. (bobspunt@gmail.com) - 2013.03.20

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, disp('USAGE: bspm_imcompare(in, option)'); return, end
if nargin<2, labels = []; end
if nargin<3, sdcut = 2.5; end

% check variable formats
if ischar(in), in = cellstr(in); end

% read in data
for i = 1:length(in)
    
    tmp = bob_reslice(in{i},in{1},1,1);
    tmp(tmp==0) = NaN;
    mim(:,i) = tmp(:);
    tmp = (tmp(:)/max(tmp(:)));
    im(:,i) = tmp(:);
    
end

if length(in)==2
    out(1) = nanmean(diff(im').^2);
    out(2) = corr(im(:,1),im(:,2),'rows','pairwise');
    fprintf('\nCOMPARING %d IMAGES\n', length(in));
    fprintf('Mean SS Difference = %2.3f\n', out(1));
    fprintf('Mean Correlation = %2.3f\n\n', out(2));
else
    M = nanmean(mim,2);
    M = (M/max(M))*100;
    out = zeros(length(in),2);
    tmpout = repmat(M,1,size(im,2)) - im;
    out(:,1) = nanmean(tmpout.^2)';
    out(:,1) = max(out(:,1)) - out(:,1); 
    for i = 1:length(in)
        out(i,2) = corr(M,im(:,i),'rows','pairwise');
    end
    zout = abs(oneoutzscore(out));
    flag = zout>sdcut;
    result = [out(:,1) flag(:,1) out(:,2) flag(:,2)];
    fprintf('\nCOMPARING %d IMAGES\n', length(in));
    disptable(result,{'MSS' 'FLAG' 'MCORR' 'FLAG'},labels,'%2.2f');
    bad = find(flag(:,1) | flag(:,2));
    disp(bad);
    disp(char(in(bad)));
    result.out = out; 
    result.idx = in(bad); 
end





% for i = 1:size(imcell,1)
%     cim = imcell(i,:);
%     csub = subname{i};
%     for p = 1:length(xpos)
%         bspm_checkreg(cim, imname, [xpos(p) 40 0], [csub ',a =  x=' num2str(xpos(p))]);
%         saveas(gcf, sprintf('%s_x=%d.pdf', csub, xpos(p)), 'pdf');
%     end
%     close all;
% end

    
 
 
 
 
