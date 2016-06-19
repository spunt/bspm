function [extent, info] = bspm_cluster_correct(im,u,alpha,range)
% BSPM_CLUSTER_CORRECT Computer extent for cluster-level correction
%
% USAGE: [k info] = bspm_cluster_correct(im,u,alpha,range)
%
%
% THIS IS A MODIFICATION OF A FUNCTION BY DRS. THOMAS NICHOLS AND MARKO
% WILKE, CorrClusTh.m. ORIGINAL DOCUMENTATION PASTED BELOW:
%
% Find the corrected cluster size threshold for a given alpha
% function [k,Pc] =CorrClusTh(SPM,u,alpha,guess)
% SPM   - SPM data structure
% u     - Cluster defining threshold
%         If less than zero, u is taken to be uncorrected P-value
% alpha - FWE-corrected level (defaults to 0.05)
% guess - Set to NaN to use a Newton-Rhapson search (default)
%         Or provide a explicit list (e.g. 1:1000) of cluster sizes to
%         search over.
%         If guess is a (non-NaN) scalar nothing happens, except the the
%         corrected P-value of guess is printed. 
%
% Finds the corrected cluster size (spatial extent) threshold for a given
% cluster defining threshold u and FWE-corrected level alpha. 
%
%_________________________________________________________________________

% -------------------------- Copyright (C) 2014 --------------------------
%	Author: Bob Spunt
%	Affilitation: Caltech
%	Email: spunt@caltech.edu
%
%	$Revision Date: Aug_20_2014

if nargin < 1
    tmp = files('*T*');
    if isempty(tmp), mfile_showhelp; return; 
    else im = tmp(1); end
end
if nargin < 2, u = .001; end
if nargin < 3, alpha = .05; end
if nargin < 4, range = 1:1000; end
if ischar(im), im = cellstr(im); end
allim = im;

for i = 1:length(allim)
    
    im = allim{i};

%% Get Design Variable %%
[impath imname] = fileparts(im);
if exist([impath filesep 'I.mat'],'file') 
    matfile = [impath filesep 'I.mat']; 
    maskfile = [impath filesep 'mask.nii'];
elseif exist([impath filesep 'SPM.mat'],'file') 
    matfile = [impath filesep 'SPM.mat'];
else
    disp('Could not find an SPM.mat or I.mat variable, exiting.'); return
end

%% Defaults %%
epsP = 1e-6;   % Corrected P-value convergence criterion (fraction of alpha)
du   = 1e-6;   % Step-size for Newton-Rhapson
maxi = 100;    % Maximum interations for refined search
STAT = 'T';    % Test statistic

%% Determime SPM or GLMFLEX %%
if strfind(matfile,'SPM.mat'), flexflag = 0; else flexflag = 1; end

%% Load and Compute Params %%
if flexflag % GLMFLEX
    load(matfile);
    try
        mask.hdr = spm_vol([I.OutputDir filesep 'mask.nii']);
    catch
        [p mf] = fileparts(im);
        mask.hdr = spm_vol([p filesep 'mask.nii']);
    end
    mask.data = spm_read_vols(mask.hdr);
    img.hdr = spm_vol(im);
    img.data = spm_read_vols(img.hdr);
    tmp = img.hdr.descrip; i1 = find(tmp=='['); i2 = find(tmp==']');
    df = str2num(tmp(i1(1)+1:i2(1)-1));
    df = [1 df];    
    n = 1;
    FWHM = I.FWHM{1};
    R = spm_resels_vol(mask.hdr,FWHM)';
    S = sum(mask.data(:)==1);
    M = I.v.mat;
    VOX  = sqrt(diag(M(1:3,1:3)'*M(1:3,1:3)))';
    FWHMmm= FWHM.*VOX; % FWHM {mm}
    v2r  = 1/prod(FWHM(~isinf(FWHM)));% voxels to resels

else % SPM
    
    load(matfile)
    df   = [1 SPM.xX.erdf];
    STAT = 'T';
    n    = 1;
    R    = SPM.xVol.R;
    S    = SPM.xVol.S;
    M    = SPM.xVol.M;
    VOX  = sqrt(diag(M(1:3,1:3)'*M(1:3,1:3)))';
    FWHM = SPM.xVol.FWHM;
    FWHMmm= FWHM.*VOX; 				% FWHM {mm}
    v2r  = 1/prod(FWHM(~isinf(FWHM))); %-voxels to resels
    
end
if ~nargout
    sf_ShowVolInfo(R,S,VOX,FWHM,FWHMmm)
end
epsP = alpha*epsP;
Status = 'OK';
if u <= 1; u = spm_u(u,df,STAT); end

if length(range)==1 & ~isnan(range)
  
  %
  % Dummy case... just report P-value
  %

  k  = range;
  Pc = spm_P(1,k*v2r,u,df,STAT,R,n,S);
  
  Status = 'JustPvalue';

elseif (spm_P(1,1*v2r,u,df,STAT,R,n,S)<alpha)

  %
  % Crazy setting, where 1 voxel cluster is significant
  %

  k = 1;
  Pc = spm_P(1,1*v2r,u,df,STAT,R,n,S);
  Status = 'TooRough';

elseif isnan(range)

  %
  % Automated search
  % 

  % Initial (lower bound) guess is the expected number of voxels per cluster
  [P Pn Em En EN] = spm_P(1,0,u,df,STAT,R,n,S);
  kr = En; % Working in resel units
  rad = (kr)^(1/3); % Parameterize proportional to cluster diameter

  %
  % Crude linear search bound answer
  %
  Pcl  = 1;   % Lower bound on P
  radu = rad; % Upper bound on rad
  Pcu  = 0;   % Upper bound on P
  radl = Inf; % Lower bound on rad
  while (Pcl > alpha)
    Pcu  = Pcl;
    radl = radu; % Save previous result
    radu = radu*1.1;
    Pcl  = spm_P(1,radu^3   ,u,df,STAT,R,n,S);
  end

  %
  % Newton-Rhapson refined search
  %
  d = 1;		    
  os = NaN;     % Old sign
  ms = (radu-radl)/10;  % Max step
  du = ms/100;
  % Linear interpolation for initial guess
  rad = radl*(alpha-Pcl)/(Pcu-Pcl)+radu*(Pcu-alpha)/(Pcu-Pcl);
  iter = 1;
  while abs(d) > epsP
    Pc  = spm_P(1,rad^3   ,u,df,STAT,R,n,S);
    Pc1 = spm_P(1,(rad+du)^3,u,df,STAT,R,n,S);
    d   = (alpha-Pc)/((Pc1-Pc)/du);
    os = sign(d);  % save old sign
    % Truncate search if step is too big
    if abs(d)>ms, 
      d = sign(d)*ms;
    end
    % Keep inside the given range
    if (rad+d)>radu, d = (radu-rad)/2; end
    if (rad+d)<radl, d = (rad-radl)/2; end
    % update
    rad = rad + d;
    iter = iter+1;
    if (iter>=maxi), 
      Status = 'TooManyIter';
      break; 
    end
  end
  % Convert back
  kr = rad^3;
  k = ceil(kr/v2r);
  Pc  = spm_P(1,k*v2r,u,df,STAT,R,n,S);

%
% Brute force!
%
else
  Pc = 1;
  for k = range
    Pc = spm_P(1,k*v2r,u,df,STAT,R,n,S);
    %fprintf('k=%d Pc=%g\n',k,Pc);
    if Pc <= alpha, 
      break; 
    end
  end;
  if (Pc > alpha)
    Status = 'OutOfRange';
  end
end

if ~nargout
    switch (Status)
     case {'JustPvalue'}
      fprintf(['  For a cluster-defining threshold of %0.4f a cluster size threshold of\n'...
           '  %d has corrected P-value %g\n\n'],...
          u,k,Pc);
     case {'OK'}
      fprintf(['  For a cluster-defining threshold of %0.4f the level %0.3f corrected\n'...
           '  cluster size threshold is %d and has size (corrected P-value) %g\n\n'],...
          u,alpha,k,Pc);
     case 'TooRough'
      fprintf(['\n  WARNING: Single voxel cluster is significant!\n\n',...
               '  For a cluster-defining threshold of %0.4f a k=1 voxel cluster\n'...
           '  size threshold has size (corrected P-value) %g\n\n'],...
          u,Pc); 
     case 'TooManyIter'
      fprintf(['\n  WARNING: Automated search failed to converge\n' ...
           '  Try systematic search.\n\n']); 
     case 'OutOfRange'  
      fprintf(['\n  WARNING: Within the range of cluster sizes searched (%g...%g)\n',...
             '  a corrected P-value <= alpha was not found (smallest P: %g)\n\n'],...
          range(1),range(end),Pc); 
      fprintf([  '  Try increasing the range or an automatic search.\n\n']); 
     otherwise
      error('Unknown status code');
    end
end
extent(i) = k;
info(i).image = im;
info(i).extent = k;
info(i).alpha = alpha;
info(i).u = u;
info(i).Pc = Pc;
end

end
function sf_ShowVolInfo(R,S,VOX,FWHM,FWHMmm)

fprintf('\n  Search Volume:  %7d %0.2fx%.02fx%0.2f voxels\n',S,VOX);
fprintf('                  %7.2f cmm, %5.2f L, %5.2f RESELS\n',...
	S*prod(VOX),S*prod(VOX)/100^3,R(end));
fprintf(['                  %0.2fx%.02fx%0.2f mm FWHM, ',...
	                  '%0.2fx%.02fx%0.2f vox FWHM\n\n'],FWHMmm,FWHM);
return
end
function P = sf_spm_get(n,wildcard,prompt)
switch spm('ver')
 case 'SPM99'
  P = spm_get(n,wildcard,prompt);
 case 'SPM2'
  if max(findstr(wildcard,'*img')) == length(wildcard)-length('*img')+1
    wildcard = [wildcard(1:end-4) '*IMAGE'];
  end
  P = spm_get(n,wildcard,prompt);
 case 'SPM5'
  wildcard = strrep(wildcard,'*','.*');
  P = spm_select(n,'mat',prompt,'',pwd,wildcard);
end
end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
