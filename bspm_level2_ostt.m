function [] = bspm_level2_ostt(cons, minN, mask, rmoutliers, omit)
% BSPM_LEVEL2_OSTT (uses GLM Flex by Aaron Schultz by default)
%
% USAGE: bspm_level2_ostt(cons, minN, mask, rmoutliers, omit)
%
%   ARGUMENTS:
%       cons: contrast images from level 1
%       name: prefix for directory in _groupstats_/<analysis_name> 
%       (optional) mask: mask file to use (default = none)
%       (optional) rmoutliers: 0 = no (default), 1 = yes
%       (optional) name: name for output directory (default is to figure this out automatically)
%       (optional) omit: string to omit subjects
%

% ------------------------------------- Copyright (C) 2014 -------------------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin<1, error('USAGE: bspm_level2_ostt(cons, minN, mask, rmoutliers, name)'); end
if nargin<2 || isempty(minN), minN = length(cons); end
if nargin<3, mask = []; end
if nargin<4, rmoutliers = 0; end
if nargin<5,
    name = sprintf('OSTT_FLEX_N=%d_minN=%d_rmout=%d', length(cons),minN, rmoutliers);
else
    if ischar(omit), omit = cellstr(omit); end
    omitidx = cellismember(cons, omit);
    if omitidx
        cons(omitidx) = [];
        omitlab = sprintf(repmat('%s_', 1, length(omitidx)), omit{:}); 
        name = sprintf('OSTT_FLEX_N=%d_minN=%d_rmout=%d_omit=%s', length(cons),minN, rmoutliers, omitlab(1:end-1));
    else
        name = sprintf('OSTT_FLEX_N=%d_minN=%d_rmout=%d', length(cons),minN, rmoutliers);
    end
end

if nargin<6, 
    namelab = {'' 'SPM_GLMFLEX'};
end

if ischar(cons), cons = cellstr(cons); end

% define output directory
if iscell(name), name = char(name); end
[p n e] = fileparts(cons{1});
idx = strfind(p,'/');
sname = p(1:idx(end-2)-1);
aname = p(idx(end)+1:length(p));
hdr = spm_vol(cons{1});
hdr = hdr.descrip;
idx = strfind(hdr,':');
cname = hdr(idx+1:end);
cname = regexprep(cname,'- All Sessions','');
cname = strtrim(cname);
gadir = fullfile(sname,'_groupstats_',aname);
gasubdir = fullfile(gadir,name);
outdir = fullfile(gasubdir,cname);
if ~isdir(gadir), mkdir(gadir); end
if ~isdir(gasubdir), mkdir(gasubdir); end
mkdir(outdir);

% setup structure
IN.N_subs = length(cons);        % Number of input images
IN.Between = 1;       % Number of levels per factor.
IN.EqualVar = 1;       % Do variance corrections
IN.Independent = 1; % Do independence correction
F = CreateDesign2(IN); 

% setup I
I.OutputDir = outdir;
I.F = F;
I.Mask = mask;
I.minN = minN;
I.DoOnlyAll = 0;
I.RemoveOutliers= rmoutliers;
I.Scans = cons;
I.Cons(1).name = cname;
I.Cons(1).Groups = 1;
I.Cons(1).Levs = 1;
I.Cons(1).ET = 1;

% run glmflex
GLM_Flex2(I);


%%% I.OutputDir = pwd;
%%% I.F = [];
%%% I.Scans = [];
%%% I.Mask = [];
%%% I.RemoveOutliers = 0;
%%% I.DoOnlyAll = 0;
%%% I.minN = 2;
%%% I.minRat = 0;
%%% I.Thresh = [];
%%% I.writeIni = 0;
%%% I.writeFin = 0;
%%% I.writeoo = 1;
%%% I.writeI = 1;
%%% I.KeepResiduals = 0;
%%% I.estSmooth = 1;
%%% I.Transform.FisherZ = 0;
%%% I.Transform.AdjustR2 = 0;
%%% I.Transform.ZScore = 0;

end

 
 
 
 
