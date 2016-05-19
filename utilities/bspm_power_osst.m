function [power,nsub] = bspm_power_osst(analysisdir, roifiles, varargin)
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

def = { ...
        'nrange',          2:100,        ...
        'alphavalue',      .05,          ...
        'desiredpower',     80 ...
        };
vals = setargs(def, varargin);
if nargin < 2, mfile_showhelp; fprintf('\t= DEFAULT SETTINGS =\n'); disp(vals); return; end
if iscell(analysisdir), analysisdir = char(analysisdir); end
d = dir(analysisdir);
c = {d(~vertcat(d.isdir)).name};
needthese = {'mask' 'ResMS' 'con_0001'};
idx2c = cellismember(c, needthese)
if ~all(idx2c)
    disp('One or more necessary files are missing:'); 
    disp(needthese(idx2c==0)); 
    return
end
ref     = fullfile(analysisdir, c(idx2c));
[d,h]   = bspm_read_vol(ref, 'reshape');
nroi    = length(roifiles);
roi     = NaN(size(d,1), nroi);
for i = 1:nroi
    tmp = bspm_reslice(roifiles(i), ref(1), 0, 1); 
    roi(:,i) = tmp(:); 
end
var = d(:,1);
con = d(:,2);
mask = d(:,3);
con(~mask) = NaN;
var(~mask) = NaN;
roi(roi==0) = NaN; 
roi(~mask,:) = NaN; 
power = zeros(length(nrange),nroi + 1);
power(:,1) = nrange; 
nsub = zeros(nroi, 1); 
for r = 1:nroi
    
    [~,roiname] = fileparts(roifiles{r}); 
    effsize = nanmean(con(roi(:,r)>0)); 
    effvar  = nanmean(sqrt(var(roi(:,r)>0)));
    str = sprintf('Effect Size = %2.4f, Variance = %2.4f', effsize, effvar); 
    printmsg(str, 'msgtitle', roiname);
    
    for n = 1:length(nrange)
        
        nsub    = nrange(n);
        c_x     = ones(nsub,1);
        c_con   = 1; 
        c_df    = nsub - 1; 
        ncp     = effsize/(effvar*sqrt(c_con*inv(c_x'*c_x)*c_con'));
        power(n,r+1)   = 100*(1-nctcdf(-1*tinv(alphavalue, c_df), c_df, ncp));
        
    end
end
if nargout==2
    for r = 1:nroi
        nsub(r) = nrange(min(find(power(:,r+1) >= desiredpower))); 
    end
end
