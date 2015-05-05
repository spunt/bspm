function [power,nsub] = bspm_power_osst(analysisdir, roifiles, nrange)
% BSPM_POWER_OSST Calculate power for one-sample t-ttest
%
%  USAGE: power = bspm_power_osst(analysisdir, roifiles, *nrange)
% __________________________________________________________________________
%  INPUTS
%	analysisdir: dir containing SPM.mat
%   roifiles: roi files
%	*nrange: subject range to calc power for (e.g., 10:100)
%

% ---------------------- Copyright (C) 2014 Bob Spunt ----------------------
%	Created:  2014-10-05
%	Email:    spunt@caltech.edu
% __________________________________________________________________________
if nargin < 2, disp('power = bspm_power_osst(analysisdir, roifiles, *nrange)'); return; end
if nargin < 3, nrange = 1:100; end
if iscell(analysisdir), analysisdir = char(analysisdir); end
needthese = {'mask.img' 'ResMS.img' 'con_0001.img'}; 
for i = 1:length(needthese)
    if ~exist(fullfile(analysisdir, needthese{i}), 'file')
        fprintf('\n%s not found in %s\n', needthese{i}, analysisdir);
        return
    end
end
alpha = .05;
desiredpower = 80; 
roi = bspm_read_vol(roifiles, 'reshape');
var = bspm_read_vol(fullfile(analysisdir, 'ResMS.img'), 'reshape'); 
con = bspm_read_vol(fullfile(analysisdir, 'con_0001.img'), 'reshape');
mask = bspm_read_vol(fullfile(analysisdir, 'mask.img'), 'reshape');
con(~mask) = NaN;
var(~mask) = NaN;
nroi = size(roi,2);
power = zeros(length(nrange),nroi);
nsub = zeros(nroi, 1); 
for r = 1:nroi
    [~,roiname] = fileparts(roifiles{r}); 
    effsize = nanmean(con(roi(:,r)>0)); 
    effvar  = nanmean(var(roi(:,r)>0));
    str = sprintf('Effect Size = %2.4f, Variance = %2.4f', effsize, effvar); 
    disp(printmsg(str, 'msgtitle', roiname))
    for nsub = nrange
        c_x     = ones(nsub,1);
        c_con   = 1; 
        c_df    = nsub - 1; 
        ncp     = effsize/(effvar*sqrt(c_con*inv(c_x'*c_x)*c_con'));
        power(nsub,r)   = 100*(1-nctcdf(-1*tinv(alpha, c_df), c_df, ncp)); 
    end
end
if nargout==2
    for r = 1:nroi
        nsub(r) = nrange(min(find(power(:,r) >= desiredpower))); 
    end
end
