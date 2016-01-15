function data = bspm_imcompare_epi_anat(epi, anat, subnames)
% BSPM_IMCOMPARE_EPI_ANAT 
%
%   USAGE: data = bspm_imcompare_epi_anat(epi, anat, subnames)
%       

% ----------------------- Copyright (C) 2014 -----------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014
if nargin<2, disp('USAGE: data = bspm_imcompare_epi_anat(epi, anat, subnames)'); return, end
if nargin<3, subnames = []; end

% check variable formats
if ischar(epi), epi = cellstr(epi); end
if ischar(anat), anat = cellstr(anat); end
if length(epi)~=length(anat), error('Different number of EPI and ANAT volumes!'); end


% read in data
for i = 1:length(epi)
    
    tmp = [];
    
    % epi
    tmpepi = bspm_reslice(epi{i},epi{1},1,1);
    tmpepi(tmpepi==0) = NaN;
    mim(:,i) = tmpepi(:);
    tmp(:,1) = (tmpepi(:)/max(tmpepi(:)));

    % anat
    tmpanat = bspm_reslice(anat{i},epi{1},1,1);
    tmpanat(tmpanat==0) = NaN;
    mimanat(:,i) = tmpanat(:);
    tmp(:,2) = (tmpanat(:)/max(tmpanat(:)));
    
    % compare
    tmp = bspm_scaledata(tmp);
    out(i,1) = nanmean(diff(tmp').^2);
    out(i,2) = corr(tmp(:,1),tmp(:,2),'rows','pairwise','type','Spearman');
    
end
zout = abs(oneoutzscore(out));
flag = zout > 2.5;
result = [out(:,1) flag(:,1) out(:,2) flag(:,2)];
fprintf('\nCOMPARING %d IMAGES\n', length(epi));
disptable(result,{'MSS' 'FLAG' 'MCORR' 'FLAG'},subnames,'%2.2f');
bad = find(flag(:,1) | flag(:,2));
disp(bad);
disp(char(epi(bad)));
data.result = result; 
data.varnames = {'MSS' 'FLAG' 'MCORR' 'FLAG'};
data.bad = bad; 
    
 
 
 
 
 
 
 
 
 
 
 
 
