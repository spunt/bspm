function bspm_level2_regression(cons, regs, regnames, mask, minN, rmoutliers)
% BSPM_LEVEL2_REGRESSION
%
%   USAGE: bspm_level2_regression(cons, regs, regname, mask, minN, rmoutliers)
%
%   ARGUMENTS:
%       cons: contrast images from level 1
%       regs: vector of values
%       regnames: names for regressors
%       mask: mask file to use (default = none)
%       minN: min number of subjects to include a voxel in the model (default = all)
%       rmoutliers: tag to remove outliers (default = 0)
%

% ------------------------------- Copyright (C) 2014 -------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<6, rmoutliers = 0; end
if nargin<5, minN = []; end
if nargin<4 || strcmp(mask,''), mask = []; maskname = 'none'; end
if nargin<3, disp('USAGE: bspm_level2_regression(cons, regs, regname, mask, minN, rmoutliers)'); return; end
if ischar(cons), cons = cellstr(cons); end
if ischar(regnames), regnames = cellstr(regnames); end
if diff(size(regnames))<0, regnames = regnames'; end
if ~isempty(mask), if iscell(mask), mask = char(mask); end; [p1 n1] = fileparts(mask); maskname = n1; end
if minN==size(regs,1), minN = []; end

%% Output Directory %%
nreg  = length(regnames);
name = sprintf(['Regress' repmat('__%s',1,nreg) '__N%dof%d_Mask=%s'],regnames{:},minN,length(cons),n1);
[p n e] = fileparts(cons{1});
idx = strfind(p,'/');
sname = p(1:idx(end-2)-1);
aname = p(idx(end)+1:length(p));
hdr = spm_vol(cons{1});
hdr = hdr.descrip;
idx = strfind(hdr,':');
cname = hdr(idx+1:end);
idx = strfind(cname,' - All');
cname = cname(1:idx-1);
cname = regexprep(cname,'- All Sessions','');
cname = strtrim(cname);
gadir = fullfile(sname,'_groupstats_',aname);
gasubdir = fullfile(gadir,name);
outdir = fullfile(gasubdir,cname);
if ~isdir(gadir), mkdir(gadir); end
if ~isdir(gasubdir), mkdir(gasubdir); end
if ~isdir(outdir), mkdir(outdir); end

%% GLM FLEX - Create Design %%
IN.Scans = cons;
IN.N_subs = [numel(cons)];
IN.Between = [1]; % constant term
IN.BetweenLabs = {{'Constant'}};
IN.EqualVar = [1]; 
IN.Independent = [1];
for r = 1:nreg
    IN.Covar{r} =  regs(:,r);
end
IN.CovarLabs = regnames;
IN.FactorLabs = [{'Const'} regnames{:}];
F = CreateDesign(IN);

%% GLM FLEX - Estimate and Compute Contrasts
I.Scans = cons;
I.OutputDir = outdir;
I.Mask = mask;
I.F = F;
if isempty(minN), I.DoOnlyAll = 1; else I.DoOnlyAll = 0; I.minN = minN; end
I.RemoveOutliers = rmoutliers;
if ~rmoutliers, I.Thresh = 0; end
I.estSmooth = 1;

%% Contrasts %%
for r = 1:nreg
    I.Cons(r).name = regnames{r};
    I.Cons(r).Groups = {1+r};
    I.Cons(r).Levs = 1;
    I.Cons(r).ET = 1;
    I.Cons(r).mean = 0;
end
if nreg > 1
    I.Cons(nreg+1).name = 'Full_Model';
    I.Cons(nreg+1).Groups = num2cell(1:nreg+1);
    I.Cons(nreg+1).Levs = 0;
    I.Cons(nreg+1).ET = 1;
    I.Cons(nreg+1).mean = 0;
end

%% Run It %%
I = GLM_Flex_Fast(I);















 
 
 
 
